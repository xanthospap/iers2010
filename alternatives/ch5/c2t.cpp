#include "iau.hpp"
#include "iersc.hpp"

iers2010::RotationMatrix3 iers2010::sofa::c2t06a(double tta, double ttb, double uta, double utb,
    double xp, double yp) noexcept {
  // Form the celestial-to-intermediate matrix for this TT.
  auto rc2i = iers2010::sofa::c2i06a(tta, ttb);

  // Predict the Earth rotation angle for this UT1.
  const double   era = iers2010::sofa::era00(uta, utb);

  // Estimate s'.
  const double   sp = iers2010::sofa::sp00(tta, ttb);

  // Form the polar motion matrix.
  auto rpom = iers2010::sofa::pom00(xp, yp, sp);

  // Combine to form the celestial-to-terrestrial matrix.
  return iers2010::sofa::c2tcio(rc2i, era, rpom);
}

iers2010::RotationMatrix3 iers2010::sofa::c2tcio(const iers2010::RotationMatrix3 &rc2i, double era, const iers2010::RotationMatrix3 &rpom) noexcept {
  auto rc2t = rc2i;
  rc2t.rotz(era);
  return rpom * rc2t;
}

iers2010::RotationMatrix3 iers2010::sofa::c2i06a(double date1, double date2) noexcept {

  // Obtain the celestial-to-true matrix (IAU 2006/2000A).
  iers2010::RotationMatrix3 rbpn;
  iers2010::sofa::pnm06a(date1, date2, rbpn);

  // Extract the X,Y coordinates.
  // iauBpn2xy(rbpn, &x, &y);
  const double x = rbpn.data[2][0];
  const double y = rbpn.data[2][1];

  // Obtain the CIO locator.
  const double s = iers2010::sofa::s06(date1, date2, x, y);

  // Form the celestial-to-intermediate matrix
  iers2010::RotationMatrix3 rc2i;
  iers2010::sofa::c2ixys(x, y, s, rc2i);

  return rc2i;
}

void iers2010::sofa::pnm06a(double date1, double date2,
            iers2010::RotationMatrix3 &rbpn) noexcept {

  double gamb, phib, psib, epsa;
  // Fukushima-Williams angles for frame bias and precession.
  iers2010::sofa::pfw06(date1, date2, gamb, phib, psib, epsa);

  double dp, de;
  // Nutation components.
  iers2010::sofa::nut06a(date1, date2, dp, de);

  // Equinox based nutation x precession x bias matrix.
  iers2010::sofa::fw2m(gamb, phib, psib + dp, epsa + de, rbpn);

  return;
}

void iers2010::sofa::fw2m(double gamb, double phib, double psi, double eps,
             iers2010::RotationMatrix3 &r) noexcept{
  r.set_identity();
  r.rotz(gamb);
  r.rotx(phib);
  r.rotz(-psi);
  r.rotx(-eps);
}

iers2010::RotationMatrix3 iers2010::sofa::pom00(double xp, double yp, double sp) noexcept {
  // initialize to identity matrix
  iers2010::RotationMatrix3 rpom;
  // apply three rotations ... W(t) = R3(−sp) x R2(xp) x R1(yp),
  rpom.rotz(sp);
  rpom.roty(-xp);
  rpom.rotx(-yp);
  return rpom;
}

double iers2010::sofa::sp00(double date1, double date2) noexcept {
  // Interval between fundamental epoch J2000.0 and current date (JC).
  double t = ((date1 - iers2010::DJ00) + date2) / iers2010::DJC;
  // Approximate s'.
  return -47e-6 * t * iers2010::DAS2R;
}
