#include "relativity.hpp"

namespace {
/* Schwarzschild terms scaled by 1/{GMe/(c^2*r^3)} */
Eigen::Matrix<double, 3, 1>
schwarzschild(const Eigen::Matrix<double, 6, 1> &sat,
              double GMe = iers2010::GMe, double beta = 1e0,
              double gamma = 1e0) noexcept {
  const auto r = sat.head<3>();
  const auto v = sat.tail<3>();
  const double s = r.norm();
  return (2e0 * (beta + gamma) * (GMe / s) - gamma * v.dot(v)) * r +
         2e0 * (1e0 + gamma) * (r.dot(v)) * v;
}

/* Lense-Thirring precession scaled by 1/{GMe/(c^2*r^3)} */
Eigen::Matrix<double, 3, 1>
lens_thirring(const Eigen::Matrix<double, 6, 1> &sat,
              [[maybe_unused]]double GMe = iers2010::GMe, double J = iers2010::Je,
              double gamma = 1e0) noexcept {
  const auto r = sat.head<3>();
  const auto v = sat.tail<3>();
  const double s = r.norm();

  /* Earth’s angular momentum per unit mass vector */
  Eigen::Matrix<double, 3, 1> j;
  j << 0e0, 0e0, J;

  return (1e0 + gamma) * ((3e0 / s / s) * r.cross(v) * r.dot(j) + v.cross(j));
}

/* Geodesic acceleration correction (de Sitter) scaled by -1/{GMe/(c^2*r^3)} */
Eigen::Matrix<double, 3, 1>
geodesic(const Eigen::Matrix<double, 6, 1> &rsat,
         const Eigen::Matrix<double, 6, 1> &rersn, double GMe = iers2010::GMe,
         double GMs = iers2010::GMs,
         double gamma = 1e0) noexcept {
  // const auto r = sat.head<3>();
  const auto v = rsat.tail<3>();
  const auto re = rersn.head<3>();
  const auto ve = rersn.tail<3>();
  const double s = rsat.head<3>().norm();
  const double R = rersn.head<3>().norm();
  return (1e0 + 2e0 * gamma) * (ve.cross(re)).cross(v) * (GMs / GMe) *
         (std::pow(s, 3) / std::pow(R, 3));
}
} /* unnamed namespace */

Eigen::Matrix<double, 3, 1> dso::iers2010_relativistic_acceleration(
    const Eigen::Matrix<double, 6, 1> &rsat,
    const Eigen::Matrix<double, 6, 1> &rsun, double GMe, double GMs, double J,
    double c, double beta, double gamma) noexcept {

  /* GMe / (c^2*r^3) */
  const double s = rsat.head<3>().norm();
  // const double f = GMe / c / c / s / s / s;
  const double f = GMe / (c * c * s * s * s);
  
  {
    auto a1 = schwarzschild(rsat, GMe, beta, gamma);
    auto a2 = lens_thirring(rsat, GMe, J, gamma);
    auto a3 = geodesic(rsat, -1e0 * rsun, GMe, GMs, gamma);
  printf("%+.3e %+.3e %+.3e\n%+.3e %+.3e %+.3e\n%+.3e %+.3e %+.3e\n%.3e\n", a1(0),a1(1),a1(2),a2(0),a2(1),a2(2),a3(0),a3(1),a3(2),f);
  }

  return (schwarzschild(rsat, GMe, beta, gamma) +
          lens_thirring(rsat, GMe, J, gamma) -
          geodesic(rsat, -1e0 * rsun, GMe, GMs, gamma)) *
         f;
}
