#include "iau.hpp"
#include "rotations.hpp"

Eigen::Matrix<double, 3, 3> iers2010::sofa::c2ixys(double x, double y,
                                                     double s) noexcept {

  /* Obtain the spherical angles E and d */
  const double r2 = x * x + y * y;
  const double e = (r2 > 0e0) ? std::atan2(y, x) : 0e0;
  const double d = std::atan(std::sqrt(r2 / (1e0 - r2)));
  
  /* R3 (−E) x R2 (−d) x R3 (E) x R3 (s),
   * X = sin d cos E,
   * Y = sin d sin E,
   * Z = cos d,
   */
  Eigen::Matrix<double, 3, 3> R = Eigen::Matrix<double, 3, 3>::Identity();
  dso::rotate<dso::RotationAxis::Z>(e,R);
  dso::rotate<dso::RotationAxis::Y>(d,R);
  dso::rotate<dso::RotationAxis::Z>(-(e+s),R);
  return R;
}
