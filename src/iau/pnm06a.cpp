#include "iau.hpp"

Eigen::Matrix<double,3,3> iers2010::sofa::pnm06a(const dso::TwoPartDate &mjd_tt) noexcept {

  double gamb, phib, psib, epsa;
  // Fukushima-Williams angles for frame bias and precession.
  iers2010::sofa::pfw06(mjd_tt, gamb, phib, psib, epsa);

  double dp, de;
  // Nutation components.
  iers2010::sofa::nut06a(mjd_tt, dp, de);

  // Equinox based nutation x precession x bias matrix.
  return iers2010::sofa::fw2m(gamb, phib, psib + dp, epsa + de);
}

