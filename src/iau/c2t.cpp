#include "iau.hpp"
#include "iersc.hpp"

#include <cstdio>

dso::Mat3x3 iers2010::sofa::c2t06a(double tta, double ttb, double uta,
                                   double utb, double xp, double yp) noexcept {

  // Form the celestial-to-intermediate matrix for this TT.
  auto rc2i = iers2010::sofa::c2i06a(tta, ttb);

  // Predict the Earth rotation angle for this UT1.
  const double era = iers2010::sofa::era00(uta, utb);

  // Estimate s'.
  const double sp = iers2010::sofa::sp00(tta, ttb);

  // Form the polar motion matrix.
  auto rpom = iers2010::sofa::pom00(xp, yp, sp);

  // Combine to form the celestial-to-terrestrial matrix.
  return iers2010::sofa::c2tcio(rc2i, era, rpom);
}

dso::Mat3x3 iers2010::sofa::c2tcio(const dso::Mat3x3 &rc2i, double era,
                                   const dso::Mat3x3 &rpom) noexcept {
  auto rc2t = rc2i;
  rc2t.rotz(era);
  return rpom * rc2t;
}

dso::Mat3x3 iers2010::sofa::c2i06a(double date1, double date2) noexcept {

  // Obtain the celestial-to-true matrix (IAU 2006/2000A).
  auto rbpn = iers2010::sofa::pnm06a(date1, date2);

  // Extract the X,Y coordinates.
  // iauBpn2xy(rbpn, &x, &y);
  const double x = rbpn.data /*[2][0]*/[6];
  const double y = rbpn.data /*[2][1]*/[7];

  // Obtain the CIO locator.
  const double s = iers2010::sofa::s06(date1, date2, x, y);

  // Form the celestial-to-intermediate matrix
  return iers2010::sofa::c2ixys(x, y, s);
}

dso::Mat3x3 iers2010::sofa::pnm06a(double date1, double date2) noexcept {

  double gamb, phib, psib, epsa;
  // Fukushima-Williams angles for frame bias and precession.
  iers2010::sofa::pfw06(date1, date2, gamb, phib, psib, epsa);

  double dp, de;
  // Nutation components.
  iers2010::sofa::nut06a(date1, date2, dp, de);

  // Equinox based nutation x precession x bias matrix.
  return iers2010::sofa::fw2m(gamb, phib, psib + dp, epsa + de);
}

dso::Mat3x3 iers2010::sofa::fw2m(double gamb, double phib, double psi,
                                 double eps) noexcept {
  dso::Mat3x3 r;
  r.rotz(gamb);
  r.rotx(phib);
  r.rotz(-psi);
  r.rotx(-eps);
  return r;
}

dso::Mat3x3 iers2010::sofa::pom00(double xp, double yp, double sp) noexcept {
  // initialize to identity matrix
  dso::Mat3x3 rpom;
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

dso::Mat3x3 iers2010::sofa::numat(double epsa, double dpsi,
                                  double deps) noexcept {
  dso::Mat3x3 rmatn;
  rmatn.rotx(epsa);
  rmatn.rotz(-dpsi);
  rmatn.rotx(-(epsa + deps));
  return rmatn;
}

dso::Mat3x3 iers2010::sofa::num06a(double date1, double date2) noexcept {
  // Mean obliquity.
  const double eps = iers2010::sofa::obl06(date1, date2);

  // Nutation components.
  double dp, de;
  iers2010::sofa::nut06a(date1, date2, dp, de);

  // Nutation matrix.
  return iers2010::sofa::numat(eps, dp, de);
}

double iers2010::sofa::eors(const dso::Mat3x3 &rnpb, double s) noexcept {
  // Evaluate Wallace & Capitaine (2006) expression (16).
  const double x = rnpb.data /*[2][0]*/[6];
  const double ax = x / (1e0 + rnpb.data /*[2][2]*/[8]);
  const double xs = 1e0 - ax * x;
  const double ys = -ax * rnpb.data /*[2][1]*/[7];
  const double zs = -x;
  const double p = rnpb.data /*[0][0]*/[0] * xs + rnpb.data /*[0][1]*/[1] * ys +
                   rnpb.data /*[0][2]*/[2] * zs;
  const double q = rnpb.data /*[1][0]*/[3] * xs + rnpb.data /*[1][1]*/[4] * ys +
                   rnpb.data /*[1][2]*/[5] * zs;
  const double eo = ((p != 0) || (q != 0)) ? s - std::atan2(q, p) : s;

  return eo;
}

double iers2010::sofa::gst06(double uta, double utb, double tta, double ttb,
                             const dso::Mat3x3 &rnpb) noexcept {
  // Extract CIP coordinates.
  // iauBpn2xy(rnpb, &x, &y);
  const double x = rnpb.data /*[2][0]*/[6];
  const double y = rnpb.data /*[2][1]*/[7];

  // The CIO locator, s.
  const double s = iers2010::sofa::s06(tta, ttb, x, y);

  // Greenwich apparent sidereal time.
  const double era = iers2010::sofa::era00(uta, utb);
  const double eors = iers2010::sofa::eors(rnpb, s);
  return iers2010::nang_02pi(era - eors);
}

void iers2010::sofa::xys00a(double date1, double date2, double &x, double &y,
                            double &s) noexcept {
  // Form the bias-precession-nutation matrix, IAU 2000A.
  auto rbpn = iers2010::sofa::pnm00a(date1, date2);

  // Extract X,Y.
  // iauBpn2xy(rbpn, x, y);
  x = rbpn.data /*[2][0]*/[6];
  y = rbpn.data /*[2][1]*/[7];

  // Obtain s.
  s = s00(date1, date2, x, y);

  return;
}

void iers2010::sofa::pr00(double date1, double date2, double &dpsipr,
                          double &depspr) noexcept {
  // Precession and obliquity corrections (radians per century)
  constexpr double PRECOR = -0.29965e0 * iers2010::DAS2R;
  constexpr double OBLCOR = -0.02524e0 * iers2010::DAS2R;

  // Interval between fundamental epoch J2000.0 and given date (JC)
  const double t = ((date1 - iers2010::DJ00) + date2) / iers2010::DJC;

  // Precession rate contributions with respect to IAU 1976/80.
  dpsipr = PRECOR * t;
  depspr = OBLCOR * t;
}

void iers2010::sofa::bi00(double &dpsibi, double &depsbi,
                          double &dra) noexcept {
  // The frame bias corrections in longitude and obliquity
  constexpr double DPBIAS = -0.041775e0 * DAS2R;
  constexpr double DEBIAS = -0.0068192e0 * DAS2R;

  // The ICRS RA of the J2000.0 equinox (Chapront et al., 2002)
  constexpr double DRA0 = -0.0146e0 * DAS2R;

  // Return the results (which are fixed).
  dpsibi = DPBIAS;
  depsbi = DEBIAS;
  dra = DRA0;
}

void iers2010::sofa::bp00(double date1, double date2, dso::Mat3x3 &rb,
                          dso::Mat3x3 &rp, dso::Mat3x3 &rbp) noexcept {
  // J2000.0 obliquity (Lieske et al. 1977)
  constexpr double EPS0 = 84'381.448e0 * iers2010::DAS2R;

  // double t, dpsibi, depsbi, dra0, psia77, oma77, chia, dpsipr, depspr, psia,
  //     oma, rbw[3][3];

  // Interval between fundamental epoch J2000.0 and current date (JC).
  const double t = ((date1 - iers2010::DJ00) + date2) / iers2010::DJC;

  // Frame bias.
  double dpsibi, depsbi, dra0;
  iers2010::sofa::bi00(dpsibi, depsbi, dra0);

  // Precession angles (Lieske et al. 1977)
  const double psia77 = (5'038.7784e0 + (-1.07259e0 + (-0.001147e0) * t) * t) *
                        t * iers2010::DAS2R;
  const double oma77 =
      EPS0 + ((0.05127e0 + (-0.007726e0) * t) * t) * t * iers2010::DAS2R;
  const double chia =
      (10.5526e0 + (-2.38064e0 + (-0.001125e0) * t) * t) * t * iers2010::DAS2R;

  // Apply IAU 2000 precession corrections.
  double dpsipr, depspr;
  iers2010::sofa::pr00(date1, date2, dpsipr, depspr);
  const double psia = psia77 + dpsipr;
  const double oma = oma77 + depspr;

  // Frame bias matrix: GCRS to J2000.0.
  rb.set_identity();
  rb.rotz(dra0);
  rb.roty(dpsibi * std::sin(EPS0));
  rb.rotx(-depsbi);

  // Precession matrix: J2000.0 to mean of date.
  rp.set_identity();
  rp.rotx(EPS0);
  rp.rotz(-psia);
  rp.rotx(-oma);
  rp.rotz(chia);

  // Bias-precession matrix: GCRS to mean of date.
  rbp = rp * rb;
}

void iers2010::sofa::pn00(double date1, double date2, double dpsi, double deps,
                          double &epsa, dso::Mat3x3 &rb, dso::Mat3x3 &rp,
                          dso::Mat3x3 &rbp, dso::Mat3x3 &rn,
                          dso::Mat3x3 &rbpn) noexcept {
  // IAU 2000 precession-rate adjustments.
  double dpsipr, depspr;
  iers2010::sofa::pr00(date1, date2, dpsipr, depspr);

  // Mean obliquity, consistent with IAU 2000 precession-nutation.
  epsa = iers2010::sofa::obl80(date1, date2) + depspr;

  // Frame bias and precession matrices and their product.
  iers2010::sofa::bp00(date1, date2, rb, rp, rbp);

  // Nutation matrix.
  rn = iers2010::sofa::numat(epsa, dpsi, deps);

  // Bias-precession-nutation matrix (classical).
  rbpn = rn * rbp;
}

void iers2010::sofa::xys06a(double date1, double date2, double &x, double &y,
                            double &s) noexcept {
  // Form the bias-precession-nutation matrix, IAU 2006/2000A
  auto rbpn = pnm06a(date1, date2);

  // Extract X,Y.
  x = rbpn.data /*[2][0]*/[6];
  y = rbpn.data /*[2][1]*/[7];
  // iauBpn2xy(rbpn, x, y);

  // Obtain s.
  s = s06(date1, date2, x, y);
}
