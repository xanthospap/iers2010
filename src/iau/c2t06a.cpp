#include "iau.hpp"

Eigen::Matrix<double, 3, 3>
iers2010::sofa::c2t06a(const dso::TwoPartDate &mjd_tt,
                       const dso::TwoPartDate &mjd_ut1, double xp,
                       double yp) noexcept {

  // Form the celestial-to-intermediate matrix for this TT.
  auto rc2i = iers2010::sofa::c2i06a(mjd_tt);

  // Predict the Earth rotation angle for this UT1.
  const double era = iers2010::sofa::era00(mjd_ut1);

  // Estimate s'.
  const double sp = iers2010::sofa::sp00(mjd_tt);

  // Form the polar motion matrix.
  auto rpom = iers2010::sofa::pom00(xp, yp, sp);

  // Combine to form the celestial-to-terrestrial matrix.
  return iers2010::sofa::c2tcio(rc2i, era, rpom);
}
