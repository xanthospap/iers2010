#include "iau.hpp"

void iers2010::sofa::xys00a(const dso::TwoPartDate &mjd_tt, double &x,
                            double &y, double &s) noexcept {
  // Form the bias-precession-nutation matrix, IAU 2000A.
  auto rbpn = iers2010::sofa::pnm00a(mjd_tt);

  // Extract X,Y.
  // iauBpn2xy(rbpn, x, y);
  x = rbpn(2, 0); //[2][0]
  y = rbpn(2, 1); //[2][1]

  // Obtain s.
  s = s00(mjd_tt, x, y);

  return;
}
