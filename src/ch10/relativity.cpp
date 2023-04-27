#include "iers2010.hpp"

namespace {
/* β PPN (parameterized post-Newtonian) parameter equal to 1 in General Relativity */
constexpr const double  rbeta = 1e0;
/* γ PPN (parameterized post-Newtonian) parameter equal to 1 in General Relativity */
constexpr const double  rgamma = 1e0;
}// unnamed namespace

/* @brief Relativistic correction to the acceleration of an artificial Earth 
 * satellite according to IERS2010 (Sec 10.3, Eq. 10.12), omitting the 
 * geodesic (de Sitter) precession correction. Acceleration in 
 * GCRS (cartesian), [m/sec]
 *
 * @param[in] rgcrs Position of the satellite with respect to the Earth in
 *                  GCRF [m]
 * @param[in] c     Speed of light [m/s ](default value is iers2010::C)
 * @param[in] J     Earth’s angular momentum per unit mass in [m^2/sec] 
 *                  (default value 9.8 × 10^8 from IERS2010)
 * @param[in] GMe   gravitational coefficient of the Earth [m^2/s^2]
 * @param[in] GMs   gravitational coefficient of the Sun [m^2/s^2]
 * @return Relativistic correction in
 */
Eigen::Matrix<double, 3, 1> iers2010::relativistic_correction(
    const Eigen::Matrix<double, 3, 1> &r, const Eigen::Matrix<double, 3, 1> &v,
    double c, double J, double gme) noexcept {
  /* norm of satellite position */
  const double s = r.norm();

  /* GMe / (c^2*r^3) */
  const double fac = gme / c / c / s / s / s;

  /* Earth’s angular momentum per unit mass vector */
  Eigen::Matrix<double, 3, 1> j;
  j << 0e0, 0e0, J;

  /* Schwarzschild terms */
  const Eigen::Matrix<double, 3, 1> Schwarzschild =
      (2e0 * (rbeta + rgamma) * (gme / s) - rgamma * v.dot(v)) * r +
      2e0 * (1e0 + rgamma) * (r.dot(v)) * v;

  /* Lense-Thirring precession */
  const Eigen::Matrix<double, 3, 1> LenseT =
      (1e0 + rgamma) * ((3e0 / s / s) * r.cross(v) * r.dot(j) + v.cross(j));

  return fac * (Schwarzschild + LenseT);
}

/* @brief Relativistic correction to the acceleration of an artificial Earth 
 * satellite according to IERS2010 (Sec 10.3, Eq. 10.12). Acceleration in 
 * GCRS (cartesian), [m/sec]
 *
 * @param[in] rgcrs Position of the satellite with respect to the Earth in
 *                  GCRF [m]
 * @param[in] Res   Position of the Earth with respect to the Sun [m]
 * @param[in] dResdt Time derivative of the position of the Earth with respect 
 *                  to the Sun [m/s]
 * @param[in] c     Speed of light [m/s ](default value is iers2010::C)
 * @param[in] J     Earth’s angular momentum per unit mass in [m^2/sec] 
 *                  (default value 9.8 × 10^8 from IERS2010)
 * @param[in] GMe   gravitational coefficient of the Earth [m^2/s^2]
 * @param[in] GMs   gravitational coefficient of the Sun [m^2/s^2]
 * @return Relativistic correction in
 */
Eigen::Matrix<double, 3, 1>
iers2010::relativistic_correction(const Eigen::Matrix<double, 3, 1> &r,
                                  const Eigen::Matrix<double, 3, 1> &Res,
                                  const Eigen::Matrix<double, 3, 1> &dResdt,
                                  const Eigen::Matrix<double, 3, 1> &v,
                                  double c, double J, double gme,
                                  double gms) noexcept {
  /* norm of satellite position */
  const double s = r.norm();

  /* GMe / (c^2*r^3) */
  const double fac = gme / c / c / s / s / s;

  /* Earth’s angular momentum per unit mass vector */
  Eigen::Matrix<double, 3, 1> j;
  j << 0e0, 0e0, J;

  /* Schwarzschild terms */
  const Eigen::Matrix<double, 3, 1> Schwarzschild =
      (2e0 * (rbeta + rgamma) * (gme / s) - rgamma * v.dot(v)) * r +
      2e0 * (1e0 + rgamma) * (r.dot(v)) * v;

  /* Lense-Thirring precession */
  const Eigen::Matrix<double, 3, 1> LenseT =
      (1e0 + rgamma) * ((3e0 / s / s) * r.cross(v) * r.dot(j) + v.cross(j));

  /* geodesic (de Sitter) precession */
  const double Resn = Res.norm();
  const auto fac2 = -gms / c / c / Resn / Resn / Resn;
  const Eigen::Matrix<double, 3, 1> DeSitter =
      (1e0 + 2 * rgamma) * ((dResdt.cross(Res)).cross(v));

  /* sum corrections */
  return fac * (Schwarzschild + LenseT) + fac2 * DeSitter;
}
