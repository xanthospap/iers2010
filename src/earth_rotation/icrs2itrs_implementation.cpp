#include "earth_rotation.hpp"
#include "iau.hpp"

namespace {

/* polar motion matrix */
auto W(double xp, double yp, double sp) noexcept {
  /* W = R3(-s') R2(xp) R1(yp) */
  return (Eigen::AngleAxisd(yp, Eigen::Vector3d::UnitX()) * 
    Eigen::AngleAxisd(xp, Eigen::Vector3d::UnitY()) * 
    Eigen::AngleAxisd(-sp, Eigen::Vector3d::UnitZ()));
}

auto C(double x, double y, double s) noexcept {
  /* Obtain the spherical angles E and d:
   * X = sin d cos E,
   * Y = sin d sin E,
   * Z = cos d
   */
  const double r2 = x * x + y * y;
  const double e = (r2 > 0e0) ? std::atan2(y, x) : 0e0;
  const double d = std::atan(std::sqrt(r2 / (1e0 - r2)));

  /* R3 (−E) x R2 (−d) x R3 (E) x R3 (s) */
  return (Eigen::AngleAxisd(e + s, Eigen::Vector3d::UnitZ()) *
          Eigen::AngleAxisd(-d, Eigen::Vector3d::UnitY()) *
          Eigen::AngleAxisd(-e, Eigen::Vector3d::UnitZ()));
}

auto R(double era) noexcept {
  return Eigen::AngleAxisd(-era, Eigen::Vector3d::UnitZ());
}

} /* unnamed namespace */

Eigen::Matrix<double, 3, 3> dso::detail::gcrs2itrs(double era, double s,
                                                   double sp, double Xcip,
                                                   double Ycip, double xp,
                                                   double yp) noexcept {
  /* compute full rotation matrix: [C x R x W] */
  auto Rm =  W(xp, yp, sp) * R(era) * C(Xcip, Ycip, s);

  /* compute rotation quaternion */
  // return dso::detail::itrs2gcrs_quaternion(
  //     era, s, sp, xcip, ycip, dso::sec2rad(meop.xp()),
  //     dso::sec2rad(meop.yp()));

  return Rm.toRotationMatrix();
}
