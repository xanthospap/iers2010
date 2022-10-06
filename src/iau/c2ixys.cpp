#include "iau.hpp"
#include "iersc.hpp"
#include <cmath>

dso::Mat3x3 iers2010::sofa::c2ixys(double x, double y, double s) noexcept {

  // Obtain the spherical angles E and d
  const double r2 = x * x + y * y;
  const double e = (r2 > 0e0) ? std::atan2(y, x) : 0e0;
  const double d = std::atan(std::sqrt(r2 / (1e0 - r2)));

  // initialize rotation matrix to identity
  dso::Mat3x3 rc2i;

  // R3 (−E) x R2 (−d) x R3 (E) x R3 (s),
  // X = sin d cos E,
  // Y = sin d sin E,
  // Z = cos d,
  rc2i.rotz(e);
  rc2i.roty(d);
  rc2i.rotz(-(e + s));

  return rc2i;
}

#ifdef USE_EIGEN
Eigen::Matrix<double, 3, 3> iers2010::sofa::c2ixys_e(double x, double y,
                                                     double s) noexcept {

  // Obtain the spherical angles E and d
  const double r2 = x * x + y * y;
  const double e = (r2 > 0e0) ? std::atan2(y, x) : 0e0;
  const double d = std::atan(std::sqrt(r2 / (1e0 - r2)));

  // initialize rotation matrix to identity
  // dso::Mat3x3 rc2i;

  // R3 (−E) x R2 (−d) x R3 (E) x R3 (s),
  // X = sin d cos E,
  // Y = sin d sin E,
  // Z = cos d,

  return Eigen::Matrix<double, 3, 3>(
      Eigen::AngleAxisd(e + s, Eigen::Vector3d::UnitZ()) *
      Eigen::AngleAxisd(-d, Eigen::Vector3d::UnitY()) *
      Eigen::AngleAxisd(-e, Eigen::Vector3d::UnitZ()));
}
#endif
