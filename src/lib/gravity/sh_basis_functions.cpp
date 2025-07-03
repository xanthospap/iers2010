#include "geodesy/transformations.hpp"
#include "gravity.hpp"
#include <array>
#include <cassert>
#include <cmath>

namespace {
/** @brief  Compute the spherical harmonic basis functions Cnm, Snm
 *
 * This function computes the spherical harmonic basis functions Cnm, Snm​ up
 * to a specified degree n, evaluated at a 3D point in space.
 *
 * These functions are the real-valued solid spherical harmonics, commonly
 * used in geodesy and gravity field modeling.
 *
 * Computes the real solid spherical harmonics (4π normalized):
 * Cnm = (1/r**(n+1)) cos(mλ) Pnm(cosθ)
 * Snm = (1/r**(n+1)) cos(mλ) Pnm(cosθ)
 *
 * At the end of the function
 * cs.C(n,m) holds the real-valued cosine basis Cnm
 * cs.S(n,m) holds the real-valued sine basis Snm​
 * These can then be multiplied by the corresponding cnm, snm gravity field
 * coefficients to compute potential, acceleration, etc.
 */
int sh_basis_cs_exterior(const Eigen::Matrix<double, 3, 1> &point,
                         dso::StokesCoeffs &cs, int max_degree,
                         int max_order) noexcept {
  /* Legendre factors:
   * Factors up to degree/order MAX_SIZE_FOR_ALF_FACTORS. Constructed only on
   * the first function call.
   */
  if (max_degree >
      dso::NormalizedLegendreFactors::MAX_SIZE_FOR_ALF_FACTORS - 3) {
    fprintf(stderr,
            "[ERROR] (Static) Size for NormalizedLegendreFactors must be "
            "augmented to perform computation (traceback: %s)\n",
            __func__);
    assert(false);
  }
  static const dso::NormalizedLegendreFactors F;

  /* clear coefficients */
  cs.clear();

  /* Normalize input point and handle scale */
  const double x = point.x() / point.squaredNorm();
  const double y = point.y() / point.squaredNorm();
  const double z = point.z() / point.squaredNorm();
  cs.C(0, 0) = 1e280 / point.norm();
  cs.S(0, 0) = 0e0;

  const double rr = point.squaredNorm();
  /* start looping, column wise */
  {
    /* column 0 */
    const int m = 0;
    /* First subdiagonal */
    cs.C(1, m) = F.f1(1, m) * z * cs.C(0, m);
    cs.S(1, m) = F.f1(1, m) * z * cs.S(0, m);
    /* rest of column */
    for (int n = m + 2; n <= max_degree; n++) {
      cs.C(n, m) =
          F.f1(n, m) * z * cs.C(n - 1, m) + F.f2(n, m) * rr * cs.C(n - 2, m);
      cs.S(n, m) =
          F.f1(n, m) * z * cs.S(n - 1, m) + F.f2(n, m) * rr * cs.S(n - 2, m);
    }
  }

  /* rest of columns, but not last column */
  for (int m = 1; m < max_order - 1; m++) {
    /* diagonal term */
    cs.C(m, m) = F.f1(m, m) * (x * cs.C(m - 1, m - 1) - y * cs.S(m - 1, m - 1));
    cs.S(m, m) = F.f1(m, m) * (y * cs.C(m - 1, m - 1) + x * cs.S(m - 1, m - 1));
    /* first sub-diagonal */
    cs.C(m + 1, m) = F.f1(m + 1, m) * z * cs.C(m, m);
    cs.S(m + 1, m) = F.f1(m + 1, m) * z * cs.S(m, m);
    /* all other rows */
    for (int n = m + 2; n <= max_degree; n++) {
      cs.C(n, m) =
          F.f1(n, m) * z * cs.C(n - 1, m) + F.f2(n, m) * rr * cs.C(n - 2, m);
      cs.S(n, m) =
          F.f1(n, m) * z * cs.S(n - 1, m) + F.f2(n, m) * rr * cs.S(n - 2, m);
    }
  }

  /* last column */
  {
    const int m = max_order;
    cs.C(m, m) = F.f1(m, m) * (x * cs.C(m - 1, m - 1) - y * cs.S(m - 1, m - 1));
    cs.S(m, m) = F.f1(m, m) * (y * cs.C(m - 1, m - 1) + x * cs.S(m - 1, m - 1));

    const int diff = max_degree - max_order;
    if (diff >= 1) {
      /* first sub-diagonal */
      cs.C(m + 1, m) = F.f1(m + 1, m) * z * cs.C(m, m);
      cs.S(m + 1, m) = F.f1(m + 1, m) * z * cs.S(m, m);
    }
    if (diff >= 2) {
      /* all other rows */
      for (int n = m + 2; n <= max_degree; n++) {
        cs.C(n, m) =
            F.f1(n, m) * z * cs.C(n - 1, m) + F.f2(n, m) * rr * cs.C(n - 2, m);
        cs.S(n, m) =
            F.f1(n, m) * z * cs.S(n - 1, m) + F.f2(n, m) * rr * cs.S(n - 2, m);
      }
    }
  }

  /* Rescale to undo 1e280 */
  cs.scale(1e-280);

  return 0;
}

/** @brief Compute Yn coefficients of the potential expansion:
 * Yn = GM/R * Σ(cnm*Cnm + snm*Snm)
 * where Cnm, Snm are the model coefficients and cnm, snm are the SH basis
 * functions.
 *
 * @param[in] point Point of computation in cartesian ECEF, [m].
 * @param[in] CS    Cnm and Snm (stokes) coefficients of the model. The degree
 * and order of the model must be less or equal to max_degree and max_order.
 * @param[out] y    The computed Yn coefficients, one per n (i.e. max_degree).
 * @param[out] cs   Scratch StokesCoeffs to perform internal computations. It
 * will be resized to hold enough elements for the computation, so it can be
 * of any size. After the call, it's size will be at least as large as
 * (max_degree, max_order), so that it can be used in subsequent calls (if
 * needed) and avoid (re)allocation.
 * @param[in] max_degree Max degree of SH expansion. max_degree should be <=
 * CS.max_degree(). If a negative number is give, then CS.max_degree() will be
 * used.
 * @param[in] max_order Max order of SH expansion. max_degree should be <=
 * CS.max_order(). If a negative number is given, then CS.max_order() will be
 * used.
 * @return Anything other than zero denotes an error.
 */
int ynm(const Eigen::Vector3d &point, const dso::StokesCoeffs &CS,
        std::vector<double> &y, dso::StokesCoeffs &cs, int max_degree = -1,
        int max_order = -1) noexcept {
  /* set/check max (n,m) */
  if (max_degree < 0)
    max_degree = CS.max_degree();
  if (max_order < 0)
    max_order = CS.max_order();

  if ((max_degree > CS.max_degree()) || (max_order >= CS.max_order())) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for Ynm computation; mismatch with "
            "given Cnm, Snm model coeffs (traceback: %s)\n",
            __func__);
    return 1;
  }

  /* resize cs (to be computed). */
  cs.resize(max_degree, max_degree);

  /* compute basis functions, cnm and snm at given point [a]*/
  if (sh_basis_cs_exterior(point / CS.Re(), cs, max_degree, max_order)) {
    fprintf(stderr, "[ERROR] Failed computing spherical harmonics basis "
                    "functions! (traceback: %s)\n");
    return 2;
  }

  /* Yn = GM/R * sum (cnm*Cnm + snm*Snm) */
  y.reserve(max_degree);
  for (int n = 0; n <= max_degree; n++) {
    double yn = cs.C(n, 0) * CS.C(n, 0);
    for (int m = 1; m <= n; m++) {
      yn += cs.C(n, m) * CS.C(n, m) + cs.S(n, m) * CS.S(n, m);
    }
    y.emplace_back(yn * (CS.GM() / CS.Re()));
  }

  return 0;
  /* Note [a], What's up with the scaling ?
   * sh_basis_cs_exterior computes cnm = (1/r^{n+1}) * Pnm * cos(ml)
   * and snm = (1/r^{n+1}) * Pnm * sin(ml)
   * but the potential reads:
   * V = (GM/r) Σ_{n} (R/n)^{n} Σ_{m} [Cnm cos(ml) Pnm + Snm sin(ml) Pnm]
   * = GM Σ_{n} R^{n} Σ_{m} [(1/r)^{n+1} * Cnm * cos(ml) * Pnm
   * + (1/r)^{n+1} * Snm * Pnm * sin(ml) ]
   * = GM Σ_{n} R^{n} Σ_{m} [Cnm * cnm + Snm * snm] (Eq. 1a)
   *
   * so, giving (point / R) in sh_basis_cs_exterior means that we will normalize
   * the radius as r' = r / R. In turn, the basis functions (cnm, snm) will be
   * computed with a scaling factor:
   *  (1/r')^{n+1} = (1/(r/R))^{n+1} = (R/r)^{n+1}.
   * In this case, Eq. 1a, becomes V = GM/R Σ_{m} [Cnm * cnm + Snm * snm]
   * which is what we actually compute (per degree n).
   */
}

int sh_deformation(const Eigen::Vector3d &point, const dso::StokesCoeffs &CS,
                   dso::StokesCoeffs &cs, Eigen::Vector3d &dr,
                   Eigen::Vector3d &gravity, double &potential,
                   int max_degree = -1, int max_order = -1) noexcept {

  /* set/check max (n,m) */
  if (max_degree < 0)
    max_degree = CS.max_degree();
  if (max_order < 0)
    max_order = CS.max_order();

  if ((max_degree > CS.max_degree()) || (max_order >= CS.max_order())) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for deformation computation; "
            "mismatch with "
            "given Cnm, Snm model coeffs (traceback: %s)\n",
            __func__);
    return 1;
  }

  /* resize cs (to be computed). */
  cs.resize(max_degree, max_degree);

  /* compute basis functions, cnm and snm at given point [a]*/
  if (sh_basis_cs_exterior(point / CS.Re(), cs, max_degree, max_order)) {
    fprintf(stderr, "[ERROR] Failed computing spherical harmonics basis "
                    "functions! (traceback: %s)\n");
    return 2;
  }

  /* displacement vector */
  dr << 0e0, 0e0, 0e0;

  /* (total) gravity at point */
  Eigen::Vector3d grav;
  grav << 0e0, 0e0, 0e0;

  /* unit vector at point */
  const Eigen::Vector3d r = point.normalized();

  for (int n = 0; n <= max_degree; n++) {
    Eigen::Vector3d am0, am1, amn;
    am0 << 0e0, 0e0, 0e0;
    am1 << 0e0, 0e0, 0e0;
    amn << 0e0, 0e0, 0e0;
    double Vm0 = 0e0, Vm1 = 0e0, Vmn = 0e0;

    /* order m=0 */
    {
      const int m = 0;

      /* potential */
      Vm0 = CS.C(n, m) * cs.C(n, m);

      /* acceleration i.e grad(V) */
      const double wm0 = std::sqrt(static_cast<double>(n + 1) * (n + 1));
      const double wp1 =
          std::sqrt(static_cast<double>(n + 1) * (n + 2)) / std::sqrt(2e0);

      const double Cm0 = wm0 * CS.C(n + 1, 0);
      const double Cp1 = wp1 * CS.C(n + 1, 1);
      const double Sp1 = wp1 * CS.S(n + 1, 1);

      const double ax = cs.C(n, 0) * (-2e0 * Cp1);
      const double ay = cs.C(n, 0) * (-2e0 * Sp1);
      const double az = cs.C(n, 0) * (-2e0 * Cm0);

      am0 = Eigen::Vector3d(ax, ay, az);
    }

    /* order m=1 */
    {
      /* only difference with the generalized formula (aka for random n,m) is in
       * wm1
       */
      const int m = 1;

      /* potential */
      Vm1 = CS.C(n, m) * cs.C(n, m) + CS.S(n, m) * cs.S(n, m);

      /* acceleration i.e grad(V) */
      const double wm1 =
          std::sqrt(static_cast<double>(n - m + 1) * (n - m + 2)) *
          std::sqrt(2e0);
      const double wm0 =
          std::sqrt(static_cast<double>(n - m + 1) * (n + m + 1));
      const double wp1 =
          std::sqrt(static_cast<double>(n + m + 1) * (n + m + 2));

      const double Cm1 = wm1 * CS.C(n + 1, m - 1);
      const double Sm1 = wm1 * CS.S(n + 1, m - 1);
      const double Cm0 = wm0 * CS.C(n + 1, m);
      const double Sm0 = wm0 * CS.S(n + 1, m);
      const double Cp1 = wp1 * CS.C(n + 1, m + 1);
      const double Sp1 = wp1 * CS.S(n + 1, m + 1);

      const double ax = cs.C(n, m) * (Cm1 - Cp1) + cs.S(n, m) * (Sm1 - Sp1);
      const double ay = cs.C(n, m) * (-Sm1 - Sp1) + cs.S(n, m) * (Cm1 + Cp1);
      const double az = cs.C(n, m) * (-2e0 * Cm0) + cs.S(n, m) * (-2e0 * Sm0);

      am1 = Eigen::Vector3d(ax, ay, az);
    }

    /* all other orders, m=2,...,max_order */
    {
      for (int m = 2; m <= std::min(n, max_order); m++) {
        /* potential */
        Vmn += CS.C(n, m) * cs.C(n, m) + CS.S(n, m) * cs.S(n, m);

        /* acceleration i.e grad(V) */
        const double wm1 =
            std::sqrt(static_cast<double>(n - m + 1) * (n - m + 2));
        const double wm0 =
            std::sqrt(static_cast<double>(n - m + 1) * (n + m + 1));
        const double wp1 =
            std::sqrt(static_cast<double>(n + m + 1) * (n + m + 2));

        const double Cm1 = wm1 * CS.C(n + 1, m - 1);
        const double Sm1 = wm1 * CS.S(n + 1, m - 1);
        const double Cm0 = wm0 * CS.C(n + 1, m);
        const double Sm0 = wm0 * CS.S(n + 1, m);
        const double Cp1 = wp1 * CS.C(n + 1, m + 1);
        const double Sp1 = wp1 * CS.S(n + 1, m + 1);

        const double ax = cs.C(n, m) * (Cm1 - Cp1) + cs.S(n, m) * (Sm1 - Sp1);
        const double ay = cs.C(n, m) * (-Sm1 - Sp1) + cs.S(n, m) * (Cm1 + Cp1);
        const double az = cs.C(n, m) * (-2 * Cm0) + cs.S(n, m) * (-2 * Sm0);

        amn += Eigen::Vector3d(ax, ay, az);
      }
    }

    const double Vn = (CS.GM() / CS.Re()) * ((Vmn + Vm1) + Vm0);
    const auto gaccn =
        CS.GM() / (2e0 * CS.Re()) *
        std::sqrt(static_cast<double>((2 * n + 1.) / (2 * n + 3.))) *
        ((amn + am1) + am0);
    grav += std::sqrt(static_cast<double>((2 * n + 1.) / (2 * n + 3.))) *
            ((amn + am1) + am0);
    dr += (hn[n] * Vn) * r + ln[n] * (gaccn - gaccn.dot(r) * r);
  } /* loop through n's */

  /* scale with gravity at point */
  grav *= CS.GM() / (2e0 * CS.Re());
  dr *= (1e0 / grav.norm());

  /* assign for output */
  gravity = grav;
  return 0;
}

} /* anonymous namespace */

template <typename C = CartesianCrd>
[[no_discard]]
int sh_basis_cs_exterior(const C &rsta, dso::StokesCoeffs &cs,
                         int max_degree = -1, int max_order = -1) noexcept {
  static_assert(dso::CoordinateTypeTraits<C>::isCartesian);

  /* set (if needed) maximum degree and order of expansion */
  if (max_degree < 0)
    max_degree = cs.max_degree();
  if (max_order < 0)
    max_order = cs.max_order();
  if ((max_order > max_degree) || (max_degree < 0 || max_order < 0)) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for spherical harmonics expansion! "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  /* check computation degree and order w.r.t. the Stokes coeffs */
  if (max_degree > cs.max_degree()) {
    fprintf(stderr,
            "[ERROR] Requesting computing SH acceleration of degree %d, but "
            "Stokes coefficients are of size %dx%d (traceback: %s)\n",
            max_degree, cs.max_degree(), cs.max_order(), __func__);
    return 1;
  }
  if (max_order > cs.max_order()) {
    fprintf(stderr,
            "[ERROR] Requesting computing SH acceleration of order %d, but "
            "Stokes coefficients are of size %dx%d (traceback: %s)\n",
            max_order, cs.max_degree(), cs.max_order(), __func__);
    return 1;
  }

  return sh_basis_cs_exterior(rsta.mv, cs, max_degree, max_order);
}

template <typename C = CartesianCrd>
[[no_discard]]
dso::StokesCoeffs sh_basis_cs_exterior(const C &rsta, int max_degree,
                                       int max_order) {
  static_assert(dso::CoordinateTypeTraits<C>::isCartesian);

  /* set (if needed) maximum degree and order of expansion */
  if ((max_order > max_degree) || (max_degree < 0 || max_order < 0)) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for spherical harmonics expansion! "
            "(traceback: %s)\n",
            __func__);
    throw std::runtime_error(err_msg);
  }

  dso::StokesCoeffs cs(max_degree, max_order);
  if (sh_basis_cs_exterior(rsta.mv, cs, max_degree, max_order)) {
    throw std::runtime_error(err_msg);
  }

  /* return the Stokes coeffs */
  return cs;
}

template <typename C = CartesianCrd>
[[no_discard]]
Eigen::Vector3d sh_deformation(const C &rsta, const dso::StokesCoeffs &CS,
                               int max_degree, int max_order,
                               Eigen::Vector3d &gravity, double &potential) {
  if (sh_deformation(rsta.mv, S, dr, gravity, potential, max_degree,
                     max_order)) {
    fprintf(stderr,
            "[ERROR] Failed computing (surface) deformation from SH "
            "expansion! (traceback: %s)\n",
            __func__);
  }
}