#include "iau.hpp"

Eigen::Matrix<double, 3, 3>
iers2010::sofa::c2i06a(const dso::TwoPartDate &mjd_tt) noexcept {

  /* Obtain the celestial-to-true matrix (IAU 2006/2000A). */
  auto rbpn = iers2010::sofa::pnm06a(mjd_tt);

  /* Extract the X,Y coordinates. */
  // iauBpn2xy(rbpn, &x, &y);
  const double x = rbpn(2, 0); // rbpn.data [2][0]
  const double y = rbpn(2, 1); // rbpn.data [2][1]

  /* Obtain the CIO locator. */
  const double s = iers2010::sofa::s06(mjd_tt, x, y);

  /* Form the celestial-to-intermediate matrix */
  return iers2010::sofa::c2ixys(x, y, s);
}
