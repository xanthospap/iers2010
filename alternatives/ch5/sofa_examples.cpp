#include "iau.hpp"
#include <cstdio>

using namespace iers2010;

/// ----------------------------------------------------------------------------
/// Examples are taken from:
/// SOFA Tools for Earth Attitude, Software version 18, Document revision 1.64
/// ----------------------------------------------------------------------------

double ex1_result[3][3] = {
  {0.973104317697512e0, 0.230363826239227e0, -0.000703163482268e0}, 
  {-0.230363800456136e0, 0.973104570632777e0, 0.000118545366806e0}, 
  {0.000711560162864e0, 0.000046626403835e0, 0.999999745754024}
};
double ex2_result[3][3] = {
  {+0.973104317697535e0, 0.230363826239128e0, -0.000703163482198},
  {-0.230363800456037e0, 0.973104570632801e0, 0.000118545366625},
  { 0.000711560162668e0, 0.000046626403995e0, 0.999999745754024}
};
double ex3_result[3][3] = {
  {0.973104317697536e0, 0.230363826239128e0, -0.000703163481769},
  {-0.230363800456036e0, 0.973104570632801e0, 0.000118545368117},
  { 0.000711560162594e0, 0.000046626402444e0, 0.999999745754024}
};
double ex4_result[3][3] = {
  {0.973104317697536e0, 0.230363826239128e0, -0.000703163481769e0},
  {-0.230363800456036e0, 0.973104570632801e0, 0.000118545368117e0},
  {0.000711560162594e0, 0.000046626402444e0, 0.999999745754024e0} 
};

void print_rotmat(const RotationMatrix3 &mat, const char *header=nullptr) noexcept {
  if (header) printf("%s\n", header);
  for (int r = 0; r < 3; r++) {
    printf("|");
    for (int j = 0; j < 3; j++) {
      printf(" %+18.15f ", mat.data[r][j]);
    }
    printf("|\n");
  }
  printf("\n");
}

RotationMatrix3 absdif(const RotationMatrix3 &mat, const double rot[3][3]) {
  RotationMatrix3 diff;
  for (int row=0; row<3; row++) {
    for (int col=0; col<3; col++) {
      diff.data[row][col] = std::abs(mat.data[row][col] - rot[row][col]);
    }
  }
  return diff;
}

RotationMatrix3 c2t06a(double tt1, double tt2, double ut11, double ut12,
                       double xp, double yp) noexcept {
  using namespace iers2010::sofa;

  // Form the celestial-to-intermediate matrix for this TT
  const auto rc2i = c2i06a(tt1, tt2);

  // Predict the Earth rotation angle for this UT1
  const double era = era00(ut11, ut12);

  // Estimate s'
  const double sp = sp00(tt1, tt2);

  // Form the polar motion matrix
  const auto rpom = pom00(xp, yp, sp);

  // combine to for the celestial-to-terrestrial matrix
  return c2tcio(rc2i, era, rpom);
}

// IAU 2000A, CIO based, using classical angles
RotationMatrix3
    iau00a_cio(double tt1, double tt2, double ut11, double ut12, double xp,
               double yp, double dx00, double dy00) noexcept {
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

  // GCRS to CIRS matrix (intermidiate result; directly computed at rc2ti)
  //auto rc2i = c2ixys(x, y, s);

  // Earth rotation angle.
  double era = era00(ut11, ut12);

  // Form celestial-terrestrial matrix (no polar motion yet).
  RotationMatrix3 rc2ti = c2ixys(x, y, s);
  rc2ti.rotz(era);
  // print_rotmat(rc2ti, "Celestial to terrestrial matrix (no polar motion)");

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
  Vector3 v1{dx00, dy00, 0e0};
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
  double gst =
      nang_02pi(gmst00(ut11, ut12, tt1, tt2) + ee00(tt1, tt2, epsa, dpsi));

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
                          double xp, double yp, double dx06,
                          double dy06) noexcept {
  // ========================= //
  // IAU 2006/2000A, CIO based //
  // ========================= //
  using namespace iers2010::sofa;

  // CIP and CIO, IAU 2006/2000A.
  double x, y, s;
  xys06a(tt1, tt2, x, y, s);

  // Add CIP corrections.
  x += dx06;
  y += dy06;

  // GCRS to CIRS matrix.
  auto rc2i = c2ixys(x, y, s);

  // Earth rotation angle.
  double era = era00(ut11, ut12);

  // Form celestial-terrestrial matrix (no polar motion yet).
  // auto rc2ti = rc2i;
  // rc2ti.rotz(era);
  rc2i.rotz(era);

  // Polar motion matrix (TIRS->ITRS, IERS 2003).
  auto sp = sp00(tt1, tt2);
  auto rpom = pom00(xp, yp, sp);

  // Form celestial-terrestrial matrix (including polar motion).
  //return rpom * rc2ti;
  return rpom * rc2i;
}

// IAU 2006/2000A, CIO based, using X,Y series
RotationMatrix3 iau06c_eq(double tt1, double tt2, double ut11, double ut12,
                          double xp, double yp, double dx06,
                          double dy06) noexcept {
  // =========================================== //
  // IAU 2006/2000A, CIO based, using X,Y series //
  // =========================================== //
  using namespace iers2010::sofa;
  
  // CIP and CIO, IAU 2006/2000A.
  double x,y;
  xy06(tt1, tt2, x, y);
  const double s = s06(tt1, tt2, x, y);
  
  // Add CIP corrections.
  x += dx06;
  y += dy06;
  
  // GCRS to CIRS matrix.
  auto rc2i = c2ixys(x, y, s);

  // Earth rotation angle.
  const double era = era00(ut11, ut12);
  
  // Form celestial-terrestrial matrix (no polar motion yet).
  rc2i.rotz(era);

  // Polar motion matrix (TIRS->ITRS, IERS 2003).
  const double sp = sp00(tt1, tt2);
  auto rpom = pom00(xp, yp, sp);
  
  // Form celestial-terrestrial matrix (including polar motion).
  return rpom * rc2i;
}

  int main() {

    // Polar motion (arcsec->radians).
    const double xp = 0.0349282 * DAS2R;
    const double yp = 0.4833163 * DAS2R;

    // UT1-UTC (s).
    // const double dut1 = -0.072073685;

    // Nutation corrections wrt IAU 1976/1980 (mas->radians).
    // const double ddp80 = -55.0655 * DMAS2R;
    // const double dde80 = -6.3580 * DMAS2R;

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

    //  IAU 2000A, CIO based, using classical angles
    auto ex1 = iau00a_cio(tt1, tt2, ut1, ut2, xp, yp, dx00, dy00);

    // IAU 2000A, equinox based, using classical angles
    auto ex2 = iau00a_eq(tt1, tt2, ut1, ut2, xp, yp, dx00, dy00);

    // IAU 2006/2000A, CIO based, using classical angles
    auto ex3 = iau06a_eq(tt1, tt2, ut1, ut2, xp, yp, dx06, dy06);
    
    // IAU 2006/2000A, CIO based, using X,Y series
    auto ex4 = iau06c_eq(tt1, tt2, ut1, ut2, xp, yp, dx06, dy06);

    printf("IAU 2000A, CIO based, using classical angles\n");
    print_rotmat(ex1);
    print_rotmat(absdif(ex1, ex1_result), "Abs. diffs from SOFA example");
    
    printf(" IAU 2000A, equinox based, using classical angles\n");
    print_rotmat(ex2);
    print_rotmat(absdif(ex2, ex2_result), "Abs. diffs from SOFA example");
    
    printf("IAU 2006/2000A, CIO based, using classical angles\n");
    print_rotmat(ex3);
    print_rotmat(absdif(ex3, ex3_result), "Abs. diffs from SOFA example");
    
    printf("IAU 2006/2000A, CIO based, using X,Y series\n");
    print_rotmat(ex4);
    print_rotmat(absdif(ex4, ex4_result), "Abs. diffs from SOFA example");

    return 0;
}