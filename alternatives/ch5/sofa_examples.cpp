#include "iau.hpp"

using namespace iers2010;

// IAU 2000A, CIO based, using classical angles
RotationMatrix3 iau00a_cio(double tt1, double tt2, double ut11, double ut12,
                           double xp, double yp, double dx00,
                           double dy00) noexcept {
  // ====================
  // IAU 2000A, CIO based
  // ====================
  using namespace iers2010::sofa;

  // CIP and CIO, IAU 2000A.
  double x, y, s;
  xys00a(tt1, tt2, x, y, s);

  // Add CIP corrections.
  x += dx00;
  y += dy00;

  // GCRS to CIRS matrix.
  auto rc2i = c2ixys(x, y, s);

  // Earth rotation angle.
  double era = era00(ut11, ut12);

  // Form celestial-terrestrial matrix (no polar motion yet).
  RotationMatrix3 rc2ti;
  rc2ti.rotz(era);

  // Polar motion matrix (TIRS->ITRS, IERS 2003).
  double sp = sp00(tt1, tt2);
  auto rpom = pom00(xp, yp, sp);

  // Form celestial-terrestrial matrix (including polar motion).
  return rpom * rc2ti;
}

// IAU 2000A, equinox based, using classical angles 
RotationMatrix3 iau00a_eq(double tt1, double tt2, double ut11, double ut12,
                           double xp, double yp, double dx00,
                           double dy00) noexcept {
  // ======================== 
  // IAU 2000A, equinox based 
  // ======================== 
  using namespace iers2010::sofa;
  
  // Nutation, IAU 2000A.
  double dp00, de00;
  nut00a(tt1, tt2, dp00, de00);
  
  // Precession-nutation quantities, IAU 2000. 
  double epsa;
  RotationMatrix3 rb, rp, rpb, rn, rnpb;
  pn00(tt1, tt2, dp00, de00, epsa, rb, rp, rpb, rn, rnpb);
  
  // Transform dX,dY corrections from GCRS to mean of date.
  Vector3 v1 {dx00, dy00, 0e0};
  Vector3 v2 = rnpb * v1;
  double ddp00 = v2.data[0] / std::sin(epsa);
  double dde00 = v2.data[1];
  
  // Corrected nutation.
  double dpsi = dp00 + ddp00;
  double deps = de00 + dde00;
  
  // Build the rotation matrix.
  rn = numat(epsa, dpsi, deps);
  
  // Combine the matrices: PN = N x P.
  // iauRxr(rn, rpb, rnpb);
  rnpb = rn * rpb;
  
  // Greenwich apparent sidereal time (IAU 2000).
  double gst = nang_02pi(gmst00(ut11, ut12, tt1, tt2) + ee00(tt1, tt2, epsa, dpsi));
  
  // Form celestial-terrestrial matrix (no polar motion yet).
  auto rc2ti = rnpb;
  rc2ti.rotz(gst);
  
  // Polar motion matrix (TIRS->ITRS, IERS 2003).
  double sp = sp00(tt1, tt2);
  auto rpom = pom00(xp, yp, sp);
  
  // Form celestial-terrestrial matrix (including polar motion).
  return rpom * rc2ti;
}

// IAU 2006/2000A, CIO based, using classical angles
RotationMatrix3 iau06a_eq(double tt1, double tt2, double ut11, double ut12,
                          double xp, double yp, double dx00,
                          double dy00) noexcept {
  // ========================= //
  // IAU 2006/2000A, CIO based //
  // ========================= //
  using namespace iers2010::sofa;

  // CIP and CIO, IAU 2006/2000A.
  double x,y,s;
  xys06a(tt1, tt2, x, y, s);
  
  // Add CIP corrections.
  x += dx06;
  y += dy06;
  
  // GCRS to CIRS matrix.
  auto rc2i = c2ixys(x, y, s);

  // Earth rotation angle.
  double era = era00(ut11, ut12);
  
  // Form celestial-terrestrial matrix (no polar motion yet).
  auto rc2ti = rc2i;
  rc2ti.rotz(era);

  // Polar motion matrix (TIRS->ITRS, IERS 2003).
  auto sp = sp00(tt1, tt2);
  auto rpom = pom00(xp, yp, sp);
  
  // Form celestial-terrestrial matrix (including polar motion).
  return rpom * rc2ti;
}


int main() {

    // Polar motion (arcsec->radians).
    const double xp = 0.0349282 * DAS2R;
    const double yp = 0.4833163 * DAS2R;

    // UT1-UTC (s).
    const double dut1 = -0.072073685;

    // Nutation corrections wrt IAU 1976/1980 (mas->radians).
    const double ddp80 = -55.0655 * DMAS2R;
    const double dde80 = -6.3580 * DMAS2R;

    // CIP offsets wrt IAU 2000A (mas->radians).
    const double dx00 = 0.1725 * DMAS2R;
    const double dy00 = -0.2650 * DMAS2R;

    // CIP offsets wrt IAU 2006/2000A (mas->radians).
    const double dx06 = 0.1750 * DMAS2R;
    const double dy06 = -0.2259 * DMAS2R;

    // date
    const double tt1 = 2454195.5e0;
    const double tt2 = 0.500754444444444e0;
    const double ut1 = tt1;
    const double ut2 = 0.499999165813831;