#include "iau.hpp"

void iers2010::sofa::xys06a(const dso::TwoPartDate &mjd_tt, double &x,
                            double &y, double &s) noexcept {
  /* Form the bias-precession-nutation matrix, IAU 2006/2000A */
  const auto rbpn = pnm06a(mjd_tt);

  /* Extract X,Y. */
  x = rbpn(2, 0); // [2][0];
  y = rbpn(2, 1); // [2][1];
  // iauBpn2xy(rbpn, x, y);

  /* Obtain s. */
  s = s06(mjd_tt, x, y);
}
