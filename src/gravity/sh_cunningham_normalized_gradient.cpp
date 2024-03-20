#include "gravity.hpp"
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fenv.h>

namespace {
/* Max size for ALF factors; if degree is more than this (-2), then it must
 * be augmented. For now, OK
 */
constexpr const int MAX_SIZE_FOR_ALF_FACTORS = 201;

struct NormalizedLegendreFactors {
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> f1;
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> f2;
  std::array<double, MAX_SIZE_FOR_ALF_FACTORS> f3;

  NormalizedLegendreFactors() noexcept
      : f1(MAX_SIZE_FOR_ALF_FACTORS, MAX_SIZE_FOR_ALF_FACTORS),
        f2(MAX_SIZE_FOR_ALF_FACTORS, MAX_SIZE_FOR_ALF_FACTORS) {
    constexpr const int N = MAX_SIZE_FOR_ALF_FACTORS;
    f1.fill_with(0e0);
    f2.fill_with(0e0);

    /* factors for the recursion ( (2n+1)/2n )^(1/2) */
    f1(1, 1) = std::sqrt(3e0);
    for (int n = 2; n < N; n++) {
      f1(n, n) = std::sqrt((2e0 * n + 1e0) / (2e0 * n));
    }

    /* factors for the recursion */
    for (int m = 0; m < N - 1; m++) {
      for (int n = m + 1; n < N; n++) {
        const double f =
            (2e0 * n + 1e0) / static_cast<double>((n + m) * (n - m));
        /* f1_nm = B_nm */
        f1(n, m) = std::sqrt(f * (2e0 * n - 1e0));
        /* f2_nm = B_nm / Bn-1m */
        f2(n, m) =
            -std::sqrt(f * (n - m - 1e0) * (n + m - 1e0) / (2e0 * n - 3e0));
      }
    }

    /* factors for acceleration */
    for (int n = 0; n < N; n++)
      f3[n] = std::sqrt((double)(2 * n + 1) / (2 * n + 3));
  }
}; /* NormalizedLegendreFactors */

/** Spherical harmonics of Earth's gravity potential to acceleration and
 *  gradient using the algorithm due to Cunningham. The acceleration and
 *  gradient are computed in Cartesian components, i.e.
 *  acceleration = (ax, ay, az), and
 *             | dax/dx dax/dy dax/dz |
 *  gradient = | day/dx day/dy day/dz |
 *             | daz/dx daz/dy daz/dz |
 *
 * @param[in] cs Normalized Stokes coefficients of spherical harmonics
 * @param[in] r  Position vector of satellite (aka point of computation) in
 *               an ECEF frame (e.g. ITRF)
 * @param[in] max_degree Max degree of spherical harmonics expansion
 * @param[in] max_order  Max order of spherical harmonics expansion
 * @param[in] Re Equatorial radius of the Earth
 * @param[in] GM Gravitational constant of Earth
 * @param[out] acc Acceleration in cartesian components in [m/s^2]
 * @param[out] gradient Gradient of acceleration in cartesian components
 * @param[in] W   A convinient storage space, as Column-Wise Lower Triangular
 *                matrix off dimensions at least (max_degree+2, max_degree+2)
 * @param[in] M   A convinient storage space, as Column-Wise Lower Triangular
 *                matrix off dimensions at least (max_degree+2, max_degree+2)
 */
int sh2gradient_cunningham_impl(
    const dso::StokesCoeffs &cs, const Eigen::Matrix<double, 3, 1> &r,
    int max_degree, int max_order, double Re, double GM,
    Eigen::Matrix<double, 3, 1> &acc, Eigen::Matrix<double, 3, 3> &gradient,
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> &W,
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise>
        &M) noexcept {

  if (max_degree > MAX_SIZE_FOR_ALF_FACTORS - 3) {
    fprintf(stderr,
            "[ERROR] (Static) Size for NormalizedLegendreFactors must be "
            "augmented to perform computation (traceback: %s)\n",
            __func__);
    assert(false);
  }

//feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  /* Factors up to degree/order MAX_SIZE_FOR_ALF_FACTORS. Constructed only on
   * the first function call
   */
  static const NormalizedLegendreFactors F;
  /* For the computations, we will need to use a larger degree (due to
   * gradient)
   */
  const int degree = max_degree;
  const int order = max_order;

  /* nullify scratch space */
  W.fill_with(0e0);
  M.fill_with(0e0);

  {
    /* (kinda) Associated Legendre Polynomials M and W
     * note that since we are going to compute derivatives (of the gravity
     * potential), we will compute M and W up to degree n+2,m+2
     */
    const Eigen::Matrix<double, 3, 1> rs = r / Re;
    const double rr = std::pow(1e0 / rs.norm(), 2);
    const double x = rs.x() * rr;
    const double y = rs.y() * rr;
    const double z = rs.z() * rr;

    /* start ALF iteration (can use scaling according to Holmes et al, 2002) */
    M(0, 0) = 1e0 / rs.norm();
    // M(0, 0) = 1e280 / rs.norm();

    /* first fill m=0 terms; note that W(n,0) = 0 (already set) */
    M(1, 0) = std::sqrt(3e0) * z * M(0, 0);
    for (int n = 2; n <= degree+2; n++) {
      M(n, 0) = F.f1(n, 0) * z * M(n - 1, 0) + F.f2(n, 0) * rr * M(n - 2, 0);
    }

    /* fill all elements for order m >= 1 */
    for (int m = 1; m < order+2; m++) {
      /* M(m,m) and W(m,m) aka, diagonal */
      M(m, m) = F.f1(m, m) * (x * M(m - 1, m - 1) - y * W(m - 1, m - 1));
      W(m, m) = F.f1(m, m) * (y * M(m - 1, m - 1) + x * W(m - 1, m - 1));

      /* if n=m+1 , we do not have a M(n-2,...) aka sub-diagonal term */
      M(m + 1, m) = F.f1(m + 1, m) * z * M(m, m);
      W(m + 1, m) = F.f1(m + 1, m) * z * W(m, m);

      /* go on .... */
      for (int n = m + 2; n <= degree+2; n++) {
        M(n, m) = F.f1(n, m) * z * M(n - 1, m) + F.f2(n, m) * rr * M(n - 2, m);
        W(n, m) = F.f1(n, m) * z * W(n - 1, m) + F.f2(n, m) * rr * W(n - 2, m);
      }
    }

    {
      /* well, we've left the lst term uncomputed */
      const int m = order+2;
      M(m, m) = F.f1(m, m) * (x * M(m - 1, m - 1) - y * W(m - 1, m - 1));
      W(m, m) = F.f1(m, m) * (y * M(m - 1, m - 1) + x * W(m - 1, m - 1));
    }

    // M.multiply(1e-280);
    // W.multiply(1e-280);

    {
      //printf("M(179) = %.15e %.15e\n", M(179,178), M(179,179));
      //printf("M(180) = %.15e %.15e %.15e\n", M(180,178), M(180,179), M(180,180));
      //printf("M(181) = %.15e %.15e %.15e %.15e\n", M(181,178), M(181,179), M(181,180), M(181,181));
      //printf("M(182) = %.15e %.15e %.15e %.15e %.15e\n", M(182,178), M(182,179), M(182,180), M(182,181), M(182,182));
      //printf("M(149) = %.15e %.15e\n", M(149,148), M(149,149));
      //printf("M(150) = %.15e %.15e %.15e\n", M(150,148), M(150,149), M(150,150));
      //printf("M(151) = %.15e %.15e %.15e %.15e\n", M(151,148), M(151,149), M(151,150), M(151,151));
      //printf("M(152) = %.15e %.15e %.15e %.15e %.15e\n", M(152,148), M(152,149), M(152,150), M(152,151), M(152,152));
    }
  } /* end computing ALF factors M and W */

  /* acceleration and gradient in cartesian components */
  acc = Eigen::Matrix<double, 3, 1>::Zero();
  gradient = Eigen::Matrix<double, 3, 3>::Zero();
  // [[maybe_unused]] const int minDegree = 0;

  /* start from smaller terms. Note that for degrees m=0,1, we are using
   * seperate loops
   */
  for (int m = order; m >= 2; --m) {
    for (int n = degree; n >= m; --n) {
      {
        /* acceleration */
        const double wm1 =
            std::sqrt(static_cast<double>(n - m + 1) * (n - m + 2));
        const double wm0 =
            std::sqrt(static_cast<double>(n - m + 1) * (n + m + 1));
        const double wp1 =
            std::sqrt(static_cast<double>(n + m + 1) * (n + m + 2));

        const double Cm1 = wm1 * M(n + 1, m - 1);
        const double Sm1 = wm1 * W(n + 1, m - 1);
        const double Cm0 = wm0 * M(n + 1, m);
        const double Sm0 = wm0 * W(n + 1, m);
        const double Cp1 = wp1 * M(n + 1, m + 1);
        const double Sp1 = wp1 * W(n + 1, m + 1);

        try {
        const double ax = cs.C(n, m) * (Cm1 - Cp1) + cs.S(n, m) * (Sm1 - Sp1);
        const double ay = cs.C(n, m) * (-Sm1 - Sp1) + cs.S(n, m) * (Cm1 + Cp1);
        const double az = cs.C(n, m) * (-2 * Cm0) + cs.S(n, m) * (-2 * Sm0);
        acc += Eigen::Matrix<double, 3, 1>(ax, ay, az) *
               std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
        } catch (std::exception &e) {
          fprintf(stderr, "FPE caught! Computation for (n,m)=(%d,%d)\n", n,m);
          fprintf(stderr, "%.15e * %.15e + %.15e * %.15e\n", cs.C(n, m), (Cm1 - Cp1), cs.S(n, m), (Sm1 - Sp1));
          fprintf(stderr, "%.15e * %.15e + %.15e * %.15e\n", cs.C(n, m), (-Sm1 - Sp1), cs.S(n, m), (Cm1 + Cp1));
          fprintf(stderr, "%.15e * %.15e + %.15e * %.15e\n", cs.C(n, m), (-2 * Cm0), cs.S(n, m), (-2 * Sm0));
          return 0;
        }

        //acc += Eigen::Matrix<double, 3, 1>(ax, ay, az) *
        //       std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
      }
      {
        /* gradient */
        const double wm2 = std::sqrt(static_cast<double>(n - m + 1) *
                                     (n - m + 2) * (n - m + 3) * (n - m + 4)) *
                           ((m == 2) ? std::sqrt(2.0) : 1.0);
        const double wm1 = std::sqrt(static_cast<double>(n - m + 1) *
                                     (n - m + 2) * (n - m + 3) * (n + m + 1));
        const double wm0 = std::sqrt(static_cast<double>(n - m + 1) *
                                     (n - m + 2) * (n + m + 1) * (n + m + 2));
        const double wp1 = std::sqrt(static_cast<double>(n - m + 1) *
                                     (n + m + 1) * (n + m + 2) * (n + m + 3));
        const double wp2 = std::sqrt(static_cast<double>(n + m + 1) *
                                     (n + m + 2) * (n + m + 3) * (n + m + 4));

        const double Cm2 = wm2 * M(n + 2, m - 2);
        const double Sm2 = wm2 * W(n + 2, m - 2);
        const double Cm1 = wm1 * M(n + 2, m - 1);
        const double Sm1 = wm1 * W(n + 2, m - 1);
        const double Cm0 = wm0 * M(n + 2, m);
        const double Sm0 = wm0 * W(n + 2, m);
        const double Cp1 = wp1 * M(n + 2, m + 1);
        const double Sp1 = wp1 * W(n + 2, m + 1);
        const double Cp2 = wp2 * M(n + 2, m + 2);
        const double Sp2 = wp2 * W(n + 2, m + 2);

        const double gxx = cs.C(n, m) * (Cm2 - 2 * Cm0 + Cp2) +
                           cs.S(n, m) * (Sm2 - 2 * Sm0 + Sp2);
        const double gxy = cs.C(n, m) * (-Sm2 + Sp2) + cs.S(n, m) * (Cm2 - Cp2);
        const double gxz = cs.C(n, m) * (-2 * Cm1 + 2 * Cp1) +
                           cs.S(n, m) * (-2 * Sm1 + 2 * Sp1);
        const double gyy = cs.C(n, m) * (-Cm2 - 2 * Cm0 - Cp2) +
                           cs.S(n, m) * (-Sm2 - 2 * Sm0 - Sp2);
        const double gyz = cs.C(n, m) * (2 * Sm1 + 2 * Sp1) +
                           cs.S(n, m) * (-2 * Cm1 - 2 * Cp1);
        const double gzz = cs.C(n, m) * (4 * Cm0) + cs.S(n, m) * (4 * Sm0);

        gradient += Eigen::Matrix<double, 3, 3>{{gxx, gxy, gxz},
                                                {gxy, gyy, gyz},
                                                {gxz, gyz, gzz}} *
                    std::sqrt((2e0 * n + 1e0) / (2e0 * n + 5e0));
      }
    } /* loop over n */
  }   /* loop over m */

  /* order m = 1 (begin summation from smaller terms) */
  for (int n = degree; n >= 1; --n) {
    const int m = 1;
    {
      /* acceleration
       * only difference with the generalized formula (aka for random n,m)
       * is in wm1
       */
      const double wm1 =
          std::sqrt(static_cast<double>(n - m + 1) * (n - m + 2)) *
          std::sqrt(2e0);
      const double wm0 =
          std::sqrt(static_cast<double>(n - m + 1) * (n + m + 1));
      const double wp1 =
          std::sqrt(static_cast<double>(n + m + 1) * (n + m + 2));

      const double Cm1 = wm1 * M(n + 1, m - 1);
      const double Sm1 = wm1 * W(n + 1, m - 1);
      const double Cm0 = wm0 * M(n + 1, m);
      const double Sm0 = wm0 * W(n + 1, m);
      const double Cp1 = wp1 * M(n + 1, m + 1);
      const double Sp1 = wp1 * W(n + 1, m + 1);

      const double ax = cs.C(n, m) * (Cm1 - Cp1) + cs.S(n, m) * (Sm1 - Sp1);
      const double ay = cs.C(n, m) * (-Sm1 - Sp1) + cs.S(n, m) * (Cm1 + Cp1);
      const double az = cs.C(n, m) * (-2e0 * Cm0) + cs.S(n, m) * (-2 * Sm0);

      acc += Eigen::Matrix<double, 3, 1>(ax, ay, az) *
             std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
    }
    {
      /* gradient */
      const double wm1 = std::sqrt(static_cast<double>(n - m + 1) *
                                   (n - m + 2) * (n - m + 3) * (n + m + 1)) *
                         std::sqrt(2e0);
      const double wm0 = std::sqrt(static_cast<double>(n - m + 1) *
                                   (n - m + 2) * (n + m + 1) * (n + m + 2));
      const double wp1 = std::sqrt(static_cast<double>(n - m + 1) *
                                   (n + m + 1) * (n + m + 2) * (n + m + 3));
      const double wp2 = std::sqrt(static_cast<double>(n + m + 1) *
                                   (n + m + 2) * (n + m + 3) * (n + m + 4));

      const double Cm1 = wm1 * M(n + 2, m - 1);
      const double Sm1 = wm1 * W(n + 2, m - 1);
      const double Cm0 = wm0 * M(n + 2, m);
      const double Sm0 = wm0 * W(n + 2, m);
      const double Cp1 = wp1 * M(n + 2, m + 1);
      const double Sp1 = wp1 * W(n + 2, m + 1);
      const double Cp2 = wp2 * M(n + 2, m + 2);
      const double Sp2 = wp2 * W(n + 2, m + 2);

      const double gxx =
          cs.C(n, m) * (-3 * Cm0 + Cp2) + cs.S(n, m) * (-Sm0 + Sp2);
      const double gxy = cs.C(n, m) * (-Sm0 + Sp2) + cs.S(n, m) * (-Cm0 - Cp2);
      const double gxz =
          cs.C(n, m) * (-2 * Cm1 + 2 * Cp1) + cs.S(n, m) * (-2 * Sm1 + 2 * Sp1);
      const double gyy =
          cs.C(n, m) * (-Cm0 - Cp2) + cs.S(n, m) * (-3 * Sm0 - Sp2);
      const double gyz =
          cs.C(n, m) * (2 * Sp1) + cs.S(n, m) * (-2 * Cm1 - 2 * Cp1);
      const double gzz = cs.C(n, m) * (4 * Cm0) + cs.S(n, m) * (4 * Sm0);

      gradient += Eigen::Matrix<double, 3, 3>{{gxx, gxy, gxz},
                                              {gxy, gyy, gyz},
                                              {gxz, gyz, gzz}} *
                  std::sqrt((2e0 * n + 1e0) / (2e0 * n + 5e0));
    }
  } /* loop over all n's for m=1 */

  /* order m = 0 */
  for (int n = degree; n >= 2; --n) {
    [[maybe_unused]] const int m = 0;
    {
      /* acceleration */
      double wm0 = std::sqrt(static_cast<double>(n + 1) * (n + 1));
      double wp1 =
          std::sqrt(static_cast<double>(n + 1) * (n + 2)) / std::sqrt(2e0);

      double Cm0 = wm0 * M(n + 1, 0);
      double Cp1 = wp1 * M(n + 1, 1);
      double Sp1 = wp1 * W(n + 1, 1);

      const double ax = cs.C(n, 0) * (-2e0 * Cp1);
      const double ay = cs.C(n, 0) * (-2e0 * Sp1);
      const double az = cs.C(n, 0) * (-2e0 * Cm0);

      acc += Eigen::Matrix<double, 3, 1>(ax, ay, az) *
             std::sqrt((2e0 * n + 1.) / (2e0 * n + 3e0));
    }
    {
      /* gradient */
      const double wm0 =
          std::sqrt(static_cast<double>(n + 1) * (n + 2) * (n + 1) * (n + 2));
      const double wp1 =
          std::sqrt(static_cast<double>(n + 1) * (n + 1) * (n + 2) * (n + 3)) /
          std::sqrt(2e0);
      const double wp2 =
          std::sqrt(static_cast<double>(n + 1) * (n + 2) * (n + 3) * (n + 4)) /
          std::sqrt(2e0);

      const double Cm0 = wm0 * M(n + 2, 0);
      const double Cp1 = wp1 * M(n + 2, 1);
      const double Sp1 = wp1 * W(n + 2, 1);
      const double Cp2 = wp2 * M(n + 2, 2);
      const double Sp2 = wp2 * W(n + 2, 2);

      const double gxx = cs.C(n, 0) * (-2e0 * Cm0 + 2e0 * Cp2);
      const double gxy = cs.C(n, 0) * (2e0 * Sp2);
      const double gxz = cs.C(n, 0) * (4e0 * Cp1);
      const double gyy = cs.C(n, 0) * (-2e0 * Cm0 - 2e0 * Cp2);
      const double gyz = cs.C(n, 0) * (4e0 * Sp1);
      const double gzz = cs.C(n, 0) * (4e0 * Cm0);

      gradient += Eigen::Matrix<double, 3, 3>{{gxx, gxy, gxz},
                                              {gxy, gyy, gyz},
                                              {gxz, gyz, gzz}} *
                  std::sqrt((2e0 * n + 1e0) / (2e0 * n + 5e0));
    }
  } /* loop over all n's for m=0 */

  /* scale ... */
  gradient *= GM / (4e0 * Re * Re * Re);
  acc *= GM / (2e0 * Re * Re);

  return 0;
}
} /* unnamed namespace */

int dso::sh2gradient_cunningham(
    const dso::StokesCoeffs &cs, const Eigen::Matrix<double, 3, 1> &r,
    Eigen::Matrix<double, 3, 1> &acc, Eigen::Matrix<double, 3, 3> &gradient,
    int max_degree, int max_order, double Re, double GM,
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> *W,
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise>
        *M) noexcept {

  /* set (if needed) maximum degree and order of expansion */
  if (max_degree < 0)
    max_degree = cs.max_degree();
  if (max_order < 0)
    max_order = cs.max_order();
  if (max_order > max_degree) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for spherical harmonics expansion! "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  /* set (if needed) geometric parameters of expansion */
  if (Re < 0)
    Re = cs.Re();
  if (GM < 0)
    GM = cs.GM();

  /* allocate (if needed) scratch space */
  int delete_mem_pool[] = {0, 0};
  if (!W) {
    W = new dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise>(
        max_degree + 2, max_degree + 2);
    delete_mem_pool[0] = 1;
  }
  if (!M) {
    M = new dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise>(
        max_degree + 2, max_degree + 2);
    delete_mem_pool[1] = 1;
  }

  /* check scratch space */
  if ((W->rows() < max_degree + 2) || (W->cols() < max_degree + 2) ||
      (M->rows() < max_degree + 2) || (M->cols() < max_degree + 2)) {
    fprintf(stderr,
            "[ERROR] Invalid size of mem pool for spherical harmonics "
            "expansion! (traceback: %s)\n",
            __func__);
    return 1;
  }

  //printf("\tCalling compute with (n,m)=(%d,%d), Re=%.9f GM=%.3f r=(%.3f %.3f "
  //       "%.3f) W(%dx%d) M(%dx%d)\n",
  //       max_degree, max_order, Re, GM, r(0), r(1), r(2), W->rows(), W->cols(), M->rows(), M->cols());

  /* call core function */
  int status = sh2gradient_cunningham_impl(cs, r, max_degree, max_order, Re, GM,
                                           acc, gradient, *W, *M);

  /* do we need to free memory ? */
  if (delete_mem_pool[0])
    delete W;
  if (delete_mem_pool[1])
    delete M;

  /* return */
  return status;
}
