#include "iau.hpp"

Eigen::Matrix<double, 3, 3> iers2010::sofa::c2ixys(double x, double y,
                                                     double s) noexcept {

  // Obtain the spherical angles E and d
  const double r2 = x * x + y * y;
  const double e = (r2 > 0e0) ? std::atan2(y, x) : 0e0;
  const double d = std::atan(std::sqrt(r2 / (1e0 - r2)));
  
  // R3 (−E) x R2 (−d) x R3 (E) x R3 (s),
  // X = sin d cos E,
  // Y = sin d sin E,
  // Z = cos d,

  return Eigen::Matrix<double, 3, 3>(
      Eigen::AngleAxisd(e + s, Eigen::Vector3d::UnitZ()) *
      Eigen::AngleAxisd(-d, Eigen::Vector3d::UnitY()) *
      Eigen::AngleAxisd(-e, Eigen::Vector3d::UnitZ()));

  /*
  { // alternative formula using equation 5.10 from IERS 2010
  const double a = .5e0 + (1e0/8e0)*(x*x + y*y);
  
  const double Q00 = 1e0 - a*x*x;
  const double Q10 = -a*x*y;
  const double Q20 = -x;

  const double Q01 = Q10;
  const double Q11 = 1e0 - a*y*y;
  const double Q21 = -y;

  const double Q02 = -Q20;
  const double Q12 = -Q21;
  const double Q22 = 1e0-a*(x*x+y*y);

  return Eigen::AngleAxisd(e, Eigen::Vector3d::UnitZ()) * Eigen::Matrix<double,3,3>(
    {Q00,Q10,Q20},{Q01,Q11,Q21},{Q02,Q12,Q22});
  }
  */
}
