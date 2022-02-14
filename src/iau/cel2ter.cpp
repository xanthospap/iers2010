#include "iau.hpp"
#include "matvec.hpp"

using dso::Mat3x3;
using dso::Vector3;

/// CGRS to ITRS via IAU 2000A, CIO based, using classical angles
/// See, SOFA doc, section 5.3
/// @param[in] tta  TT as a 2-part Julian Date. The TT and UT1 dates tta+ttb
///            and uta+utb are Julian Dates, apportioned in any convenient way
///            between the arguments uta and utb.
/// @param[in] ttb  TT as a 2-part Julian Date (see above)
/// @param[in] uta  UT1 as a 2-part Julian Date. In the case of uta,utb, the
///            date & time method is best matched to the Earth rotation angle
///            algorithm used:  maximum precision is delivered when the uta
///            argument is for 0hrs UT1 on the day in question and the utb
///            argument lies in the range 0 to 1, or vice versa.
/// @param[in] utb  UT1 as a 2-part Julian Date (see above)
/// @param[in] xp X-coordinate of the pole (radians). The arguments xp and yp
///            are the coordinates (in radians) of the Celestial Intermediate
///            Pole with respect to the International Terrestrial Reference
///            System (see IERS Conventions 2003), measured along the meridians
///            0 and 90 deg west respectively.
/// @param[in] yp Y-coordinate of the pole (radians). See above
/// @param[in] dx00 and dy00 are corrections to the celestial pole with respect
///            to the IAU 2000A model (published by IERS or in ERP noted as
///            dX_2000 and dY_2000). Units are radians.
/// @param[in] dy00 corrections to the celestial pole, see above
Mat3x3 gcrs2itrs_00aCio(double tta, double ttb, double ut1a,
                                           double ut1b, double xp, double yp,
                                           double dx00, double dy00) noexcept {
  // ====================
  // IAU 2000A, CIO based
  // ====================
  using namespace iers2010::sofa;

  // CIP and CIO, IAU 2000A.
  double x, y, s;
  xys00a(tta, ttb, x, y, s);

  // Add CIP corrections.
  x += dx00;
  y += dy00;

  // GCRS to CIRS matrix (intermidiate result; directly computed at rc2ti)
  // auto rc2i = c2ixys(x, y, s);

  // Earth rotation angle.
  const double era = era00(ut1a, ut1b);

  // Form celestial-terrestrial matrix (no polar motion yet).
  Mat3x3 rc2ti = c2ixys(x, y, s);
  rc2ti.rotz(era);

  // Polar motion matrix (TIRS->ITRS, IERS 2003).
  const double sp = sp00(tta, ttb);
  const auto rpom = pom00(xp, yp, sp);

  // Form celestial-terrestrial matrix (including polar motion).
  return rpom * rc2ti;
}

/// CGRS to ITRS via IAU 2000A, equinox based, using classical angles
/// See, SOFA doc, section 5.4
/// @param[in] tta  TT as a 2-part Julian Date. The TT and UT1 dates tta+ttb
///            and uta+utb are Julian Dates, apportioned in any convenient way
///            between the arguments uta and utb.
/// @param[in] ttb  TT as a 2-part Julian Date (see above)
/// @param[in] uta  UT1 as a 2-part Julian Date. In the case of uta,utb, the
///            date & time method is best matched to the Earth rotation angle
///            algorithm used:  maximum precision is delivered when the uta
///            argument is for 0hrs UT1 on the day in question and the utb
///            argument lies in the range 0 to 1, or vice versa.
/// @param[in] utb  UT1 as a 2-part Julian Date (see above)
/// @param[in] xp X-coordinate of the pole (radians). The arguments xp and yp
///            are the coordinates (in radians) of the Celestial Intermediate
///            Pole with respect to the International Terrestrial Reference
///            System (see IERS Conventions 2003), measured along the meridians
///            0 and 90 deg west respectively.
/// @param[in] yp Y-coordinate of the pole (radians). See above
/// @param[in] dx00 and dy00 are corrections to the celestial pole with respect
///            to the IAU 2000A model (published by IERS or in ERP noted as
///            dX_2000 and dY_2000). Units are radians.
/// @param[in] dy00 corrections to the celestial pole, see above
Mat3x3
gcrs2itrs_00a_equinox(double tta, double ttb, double ut1a, double ut1b,
                      double xp, double yp, double dx00, double dy00) noexcept {
  // ========================
  // IAU 2000A, equinox based
  // ========================
  using namespace iers2010::sofa;

  // Nutation, IAU 2000A.
  double dp00, de00;
  nut00a(tta, ttb, dp00, de00);

  // Precession-nutation quantities, IAU 2000.
  double epsa;
  Mat3x3 rb, rp, rpb, rn, rnpb;
  pn00(tta, ttb, dp00, de00, epsa, rb, rp, rpb, rn, rnpb);

  // Transform dX,dY corrections from GCRS to mean of date.
  Vector3 v1{dx00, dy00, 0e0};
  Vector3 v2 = rnpb * v1;
  const double ddp00 = v2.data[0] / std::sin(epsa);
  const double dde00 = v2.data[1];

  // Corrected nutation.
  const double dpsi = dp00 + ddp00;
  const double deps = de00 + dde00;

  // Build the rotation matrix.
  rn = numat(epsa, dpsi, deps);

  // Combine the matrices: PN = N x P.
  // iauRxr(rn, rpb, rnpb);
  rnpb = rn * rpb;

  // Greenwich apparent sidereal time (IAU 2000).
  const double gst = iers2010::nang_02pi(gmst00(ut1a, ut1b, tta, ttb) +
                                   ee00(tta, ttb, epsa, dpsi));

  // Form celestial-terrestrial matrix (no polar motion yet).
  auto rc2ti = rnpb;
  rc2ti.rotz(gst);

  // Polar motion matrix (TIRS->ITRS, IERS 2003).
  const double sp = sp00(tta, ttb);
  const auto rpom = pom00(xp, yp, sp);

  // Form celestial-terrestrial matrix (including polar motion).
  return rpom * rc2ti;
}

/// CGRS to ITRS via IAU 2006/2000A, CIO based, using classical angles
/// See, SOFA doc, section 5.5
/// @param[in] tta  TT as a 2-part Julian Date. The TT and UT1 dates tta+ttb
///            and uta+utb are Julian Dates, apportioned in any convenient way
///            between the arguments uta and utb.
/// @param[in] ttb  TT as a 2-part Julian Date (see above)
/// @param[in] uta  UT1 as a 2-part Julian Date. In the case of uta,utb, the
///            date & time method is best matched to the Earth rotation angle
///            algorithm used:  maximum precision is delivered when the uta
///            argument is for 0hrs UT1 on the day in question and the utb
///            argument lies in the range 0 to 1, or vice versa.
/// @param[in] utb  UT1 as a 2-part Julian Date (see above)
/// @param[in] xp X-coordinate of the pole (radians). The arguments xp and yp
///            are the coordinates (in radians) of the Celestial Intermediate
///            Pole with respect to the International Terrestrial Reference
///            System (see IERS Conventions 2003), measured along the meridians
///            0 and 90 deg west respectively.
/// @param[in] yp Y-coordinate of the pole (radians). See above
/// @param[in] dx06 and dy00 are corrections to the celestial pole with respect
///            to the IAU 2006/2000A model. Units are radians.
/// @param[in] dy06 corrections to the celestial pole, see above
Mat3x3 gcrs2itrs_06_cio(double tta, double ttb, double ut1a,
                                           double ut1b, double xp, double yp,
                                           double dx06, double dy06) noexcept {
  // ========================= //
  // IAU 2006/2000A, CIO based //
  // ========================= //
  using namespace iers2010::sofa;

  // CIP and CIO, IAU 2006/2000A.
  double x, y, s;
  xys06a(tta, ttb, x, y, s);

  // Add CIP corrections.
  x += dx06;
  y += dy06;

  // Earth rotation angle.
  const double era = era00(ut1a, ut1b);
  
  // GCRS to CIRS matrix.
  auto rc2i = c2ixys(x, y, s);

  // Form celestial-terrestrial matrix (no polar motion yet).
  // auto rc2ti = rc2i;
  // rc2ti.rotz(era);
  rc2i.rotz(era);

  // Polar motion matrix (TIRS->ITRS, IERS 2003).
  const auto sp = sp00(tta, ttb);
  const auto rpom = pom00(xp, yp, sp);

  // Form celestial-terrestrial matrix (including polar motion).
  // return rpom * rc2ti;
  return rpom * rc2i;
}

/// CGRS to ITRS via IAU 2006/2000A, CIO based, using X, Y series
/// See, SOFA doc, section 5.6
/// @param[in] tta  TT as a 2-part Julian Date. The TT and UT1 dates tta+ttb
///            and uta+utb are Julian Dates, apportioned in any convenient way
///            between the arguments uta and utb.
/// @param[in] ttb  TT as a 2-part Julian Date (see above)
/// @param[in] uta  UT1 as a 2-part Julian Date. In the case of uta,utb, the
///            date & time method is best matched to the Earth rotation angle
///            algorithm used:  maximum precision is delivered when the uta
///            argument is for 0hrs UT1 on the day in question and the utb
///            argument lies in the range 0 to 1, or vice versa.
/// @param[in] utb  UT1 as a 2-part Julian Date (see above)
/// @param[in] xp X-coordinate of the pole (radians). The arguments xp and yp
///            are the coordinates (in radians) of the Celestial Intermediate
///            Pole with respect to the International Terrestrial Reference
///            System (see IERS Conventions 2003), measured along the meridians
///            0 and 90 deg west respectively.
/// @param[in] yp Y-coordinate of the pole (radians). See above
/// @param[in] dx06 and dy00 are corrections to the celestial pole with respect
///            to the IAU 2006/2000A model. Units are radians.
/// @param[in] dy06 corrections to the celestial pole, see above
Mat3x3 gcrs2itrs_06_xys(double tta, double ttb, double ut1a,
                                           double ut1b, double xp, double yp,
                                           double dx06, double dy06) noexcept {
  // =========================================== //
  // IAU 2006/2000A, CIO based, using X,Y series //
  // =========================================== //
  using namespace iers2010::sofa;

  // CIP and CIO, IAU 2006/2000A.
  double x, y;
  xy06(tta, ttb, x, y);
  const double s = s06(tta, ttb, x, y);

  // Add CIP corrections.
  x += dx06;
  y += dy06;

  // Earth rotation angle.
  const double era = era00(ut1a, ut1b);

  // GCRS to CIRS matrix.
  auto rc2i = c2ixys(x, y, s);

  // Form celestial-terrestrial matrix (no polar motion yet).
  rc2i.rotz(era);

  // Polar motion matrix (TIRS->ITRS, IERS 2003).
  const double sp = sp00(tta, ttb);
  auto rpom = pom00(xp, yp, sp);

  // Form celestial-terrestrial matrix (including polar motion).
  return rpom * rc2i;
}
