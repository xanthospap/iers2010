#include "iau.hpp"

Eigen::Matrix<double, 3, 3>
iers2010::sofa::num06a(const dso::TwoPartDate &mjd_tt) noexcept {
  /* Mean obliquity. */
  const double eps = iers2010::sofa::obl06(mjd_tt);

  /* Nutation components */
  double dp, de;
  iers2010::sofa::nut06a(mjd_tt, dp, de);

  /* Nutation matrix */
  return iers2010::sofa::numat(eps, dp, de);
}
