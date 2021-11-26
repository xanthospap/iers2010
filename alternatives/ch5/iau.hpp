#ifndef __IAU_IERS10_SOFA_CPP_HPP__
#define __IAU_IERS10_SOFA_CPP_HPP__

#include "iersc.hpp"
#include <cmath>
#include <cstring>

/// @note Note on 2-part dates:
/// The whatever date dj1+dj2 is a Julian Date, apportioned in any
//  convenient way between the arguments dj1 and dj2.  For example,
//  JD(UT1)=2450123.7 could be expressed in any of these ways,
//  among others:
//
//            dj1            dj2
//
//        2450123.7           0.0       (JD method)
//        2451545.0       -1421.3       (J2000 method)
//        2400000.5       50123.2       (MJD method)
//        2450123.5           0.2       (date & time method)

namespace iers2010 {

struct RotationMatrix3 {
  double data[3][3] = {{1e0, 0e0, 0e0}, {0e0, 1e0, 0e0}, {0e0, 0e0, 1e0}};
  
  RotationMatrix3(){};
  void memcp(const RotationMatrix3 &a) noexcept {
    std::memcpy(data[0], a.data[0], sizeof(double) * 3);
    std::memcpy(data[1], a.data[1], sizeof(double) * 3);
    std::memcpy(data[2], a.data[2], sizeof(double) * 3);
  }
  
  RotationMatrix3 &operator=(const RotationMatrix3 &a) noexcept {
    if (this != &a)
      this->memcp(a);
    return *this;
  }
  
  RotationMatrix3 &operator=(RotationMatrix3 &&a) noexcept {
    this->memcp(a);
    return *this;
  }
  
  RotationMatrix3(const RotationMatrix3 &a) noexcept { this->memcp(a); }
  
  RotationMatrix3(RotationMatrix3 &&a) noexcept { this->memcp(a); }
  
  /// @brief Multiply two 3x3 matrices (aka this * b)
  RotationMatrix3 operator*(const RotationMatrix3 &b) const noexcept;
  
  /// @brief Multiply two 3x3 matrices (aka this * b) and store result in
  ///        this instance
  void mult_inplace(const RotationMatrix3 &b) noexcept;
  
  /// @brief Transpose a 3x3 matric (in place)
  void transpose_inplace() noexcept;
  
  /// @brief Set to identity matrix
  void set_identity() noexcept;
  
  /// @brief Rotate an r-matrix about the x-axis.
  /// This will actually perform the operation R = Rx * R, with Rx =
  ///  1        0            0
  ///  0   + cos(phi)   + sin(phi)
  ///  0   - sin(phi)   + cos(phi)
  /// @param[in] angle (radians)
  void rotx(double angle) noexcept;
  
  /// @brief Rotate an r-matrix about the y-axis.
  /// This will actually perform the operation R = Ry * R, with Rx =
  ///  + cos(phi)     0      - sin(phi)
  ///       0           1           0
  ///  + sin(phi)     0      + cos(phi)
  /// @param[in] angle (radians)
  void roty(double) noexcept;
  
  /// @brief Rotate an r-matrix about the z-axis.
  /// This will actually perform the operation R = Rz * R, with Rx =
  ///  + cos(psi)   + sin(psi)     0
  ///  - sin(psi)   + cos(psi)     0
  ///       0            0         1
  /// @param[in] angle (radians)
  void rotz(double) noexcept;
}; // RotationMatrix3

namespace sofa {

/// @brief Formulate celestial to terrestrial matrix 
/// Form the celestial to terrestrial matrix given the date, the UT1 and the 
/// polar motion, using the IAU 2006/2000A precession-nutation model.
/// The matrix rc2t transforms from celestial to terrestrial
/// coordinates:
///
///        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
///
///              = rc2t * [CRS]
///
/// where [CRS] is a vector in the Geocentric Celestial Reference System and 
/// [TRS] is a vector in the International Terrestrial Reference System (see 
/// IERS Conventions 2003), RC2I is the celestial-to-intermediate matrix, ERA 
/// is the Earth rotation angle and RPOM is the polar motion matrix.
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
/// @param[in] xp   X-coordinate of the pole. The arguments xp and yp are the 
///            coordinates (in radians) of the Celestial Intermediate Pole with 
///            respect to the International Terrestrial Reference System (see 
///            IERS Conventions 2003), measured along the meridians 0 and 90 deg 
///            west respectively.
/// @param[in] yp coordinates of the pole (radians, see above)
RotationMatrix3 c2t06a(double tta, double ttb, double uta, double utb,
               double xp, double yp) noexcept;

/// @brief Formulate celestial to terrestrial matrix
/// Assemble the celestial to terrestrial matrix from CIO-based components (the 
/// celestial-to-intermediate matrix, the Earth Rotation Angle and the polar 
/// motion matrix).
/// The relationship between the arguments is as follows:
///
///       [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
///             = rc2t * [CRS]
///
/// where [CRS] is a vector in the Geocentric Celestial Reference System and 
/// [TRS] is a vector in the International Terrestrial Reference System (see 
/// IERS Conventions 2003).
/// @param[in] rc2i celestial-to-intermediate matrix
/// @param[in] era  Earth rotation angle (radians)
/// @param[in] rpom polar-motion matrix
/// @return rc2t celestial-to-terrestrial matrix
RotationMatrix3 c2tcio(const RotationMatrix3 &rc2i, double era, const RotationMatrix3 &rpom) noexcept;

/// @brief Celestial-to-intermediate matrix
/// Form the celestial-to-intermediate matrix for a given date using the IAU 
/// 2006 precession and IAU 2000A nutation models.
/// The matrix rc2i is the first stage in the transformation from celestial to 
/// terrestrial coordinates:
/// [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]
///        =  RC2T * [CRS]
/// where [CRS] is a vector in the Geocentric Celestial Reference System and 
/// [TRS] is a vector in the International Terrestrial Reference System (see 
/// IERS Conventions 2003), ERA is the Earth Rotation Angle and RPOM is the 
/// polar motion matrix.
/// @param[in] date1 (date2) TT as a 2-part Julian Date. The TT date 
///            date1+date2 is a Julian Date, apportioned in any convenient way 
///            between the two arguments. Optimally, the 'J2000 method' is best 
///            matched to the way the argument is handled internally and will 
///            deliver the optimum resolution.
/// @param[in] date2 (date1) TT as a 2-part Julian Date. See above
RotationMatrix3 c2i06a(double date1, double date2) noexcept;

/// @brief Form precession-nutation matrix
/// Form the matrix of precession-nutation for a given date (including frame 
/// bias), equinox based, IAU 2006 precession and IAU 2000A nutation models.
/// The matrix operates in the sense V(date) = rbpn * V(GCRS), where the 
/// p-vector V(date) is with respect to the true equatorial triad of date 
/// date1+date2 and the p-vector V(GCRS) is with respect to the Geocentric 
/// Celestial Reference System (IAU, 2000).
/// @param[in] date1 (date2) TT as a 2-part Julian Date. The TT date 
///            date1+date2 is a Julian Date, apportioned in any convenient way 
///            between the two arguments. Optimally, the 'J2000 method' is best 
///            matched to the way the argument is handled internally and will 
///            deliver the optimum resolution.
/// @param[in] date2 (date1) TT as a 2-part Julian Date. See above
/// @param[out] rbpn bias-precession-nutation matrix
void pnm06a(double date1, double date2,
            iers2010::RotationMatrix3 &rbpn) noexcept;

/// @brief Form rotation matrix given the Fukushima-Williams angles.
/// 1) Naming the following points:
///
///           e = J2000.0 ecliptic pole,
///           p = GCRS pole,
///           E = ecliptic pole of date,
///     and   P = CIP,
///
///     the four Fukushima-Williams angles are as follows:
///
///        gamb = gamma = epE
///        phib = phi = pE
///        psi = psi = pEP
///        eps = epsilon = EP
///
///  2) The matrix representing the combined effects of frame bias,
///     precession and nutation is:
///
///        NxPxB = R_1(-eps).R_3(-psi).R_1(phib).R_3(gamb)
///
///  3) The present function can construct three different matrices,
///     depending on which angles are supplied as the arguments gamb,
///     phib, psi and eps:
///
///     o  To obtain the nutation x precession x frame bias matrix,
///        first generate the four precession angles known conventionally
///        as gamma_bar, phi_bar, psi_bar and epsilon_A, then generate
///        the nutation components Dpsi and Depsilon and add them to
///        psi_bar and epsilon_A, and finally call the present function
///        using those four angles as arguments.
///
///     o  To obtain the precession x frame bias matrix, generate the
///        four precession angles and call the present function.
///
///     o  To obtain the frame bias matrix, generate the four precession
///        angles for date J2000.0 and call the present function.
///
///     The nutation-only and precession-only matrices can if necessary
///     be obtained by combining these three appropriately.
/// @param[in]   gamb F-W angle gamma_bar (radians)
/// @param[in]   phib F-W angle phi_bar (radians)
/// @param[in]   psi  F-W angle psi (radians)
/// @param[in]   eps  F-W angle epsilon (radians)
void fw2m(double gamb, double phib, double psi, double eps,
             iers2010::RotationMatrix3 &r) noexcept;

/// @brief IAU 2000A nutation with adjustments to match the IAU 2006 precession.
/// The nutation components in longitude and obliquity are in radians and with 
/// respect to the mean equinox and ecliptic of date, IAU 2006 precession model 
/// (Hilton et al. 2006, Capitaine et al. 2005).
/// @param[in] date1 (date2) TT as a 2-part Julian Date. The TT date 
///            date1+date2 is a Julian Date, apportioned in any convenient way 
///            between the two arguments. Optimally, the 'J2000 method' is best 
///            matched to the way the argument is handled internally and will 
///            deliver the optimum resolution.
/// @param[in] date2 (date1) TT as a 2-part Julian Date. See above
/// @param[out] dpsi nutation (luni-solar+planetary) longtitude component in 
///             radians
/// @param[out] deps nutation (luni-solar+planetary) obliquity component in 
///             radians
void nut06a(double date1, double date2, double &dpsi, double &deps) noexcept;

/// @brief Nutation, IAU 2000A model
/// Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation with 
/// free core nutation omitted).
/// The nutation components in longitude and obliquity are in radians and with 
/// respect to the equinox and ecliptic of date.  The obliquity at J2000.0 is 
/// assumed to be the Lieske et al. (1977) value of 84381.448 arcsec.
/// 
/// Both the luni-solar and planetary nutations are included. The latter are 
/// due to direct planetary nutations and the perturbations of the lunar and 
/// terrestrial orbits.
///
/// The function computes the MHB2000 nutation series with the associated 
/// corrections for planetary nutations. It is an implementation of the nutation 
/// part of the IAU 2000A precession-nutation model, formally adopted by the 
/// IAU General Assembly in 2000, namely MHB2000 (Mathews et al. 2002), but 
/// with the free core nutation (FCN - see Note below) omitted.
///
/// The full MHB2000 model also contains contributions to the nutations in 
/// longitude and obliquity due to the free-excitation of the free-core-nutation 
/// during the period 1979-2000.  These FCN terms, which are time-dependent and 
/// unpredictable, are NOT included in the present function and, if required, 
/// must be independently computed.  With the FCN corrections included, the
/// present function delivers a pole which is at current epochs accurate to a 
/// few hundred microarcseconds.  The omission of FCN introduces further errors 
/// of about that size.
/// 
/// The present function provides classical nutation.  The MHB2000 algorithm, 
/// from which it is adapted, deals also with (i) the offsets between the GCRS 
/// and mean poles and (ii) the adjustments in longitude and obliquity due to 
/// the changed precession rates. These additional functions, namely frame bias 
/// and precession adjustments, are supported by the SOFA functions bi00 
/// (iauBi00) and pr00 (iauPr00).
///
/// The MHB2000 algorithm also provides "total" nutations, comprising the 
/// arithmetic sum of the frame bias, precession adjustments, luni-solar 
/// nutation and planetary nutation.  These total nutations can be used in 
/// combination with an existing IAU 1976 precession implementation, such as 
/// iauPmat76,  to deliver GCRS-to-true predictions of sub-mas accuracy at 
/// current dates. However, there are three shortcomings in the MHB2000 model 
/// that must be taken into account if more accurate or definitive results
/// are required (see Wallace 2002):
///
///       (i) The MHB2000 total nutations are simply arithmetic sums,
///           yet in reality the various components are successive Euler
///           rotations.  This slight lack of rigor leads to cross terms
///           that exceed 1 mas after a century.  The rigorous procedure
///           is to form the GCRS-to-true rotation matrix by applying the
///           bias, precession and nutation in that order.
///
///      (ii) Although the precession adjustments are stated to be with
///           respect to Lieske et al. (1977), the MHB2000 model does
///           not specify which set of Euler angles are to be used and
///           how the adjustments are to be applied.  The most literal
///           and straightforward procedure is to adopt the 4-rotation
///           epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR
///           to psi_A and DEPSPR to both omega_A and eps_A.
///
///     (iii) The MHB2000 model predates the determination by Chapront
///           et al. (2002) of a 14.6 mas displacement between the
///           J2000.0 mean equinox and the origin of the ICRS frame.  It
///           should, however, be noted that neglecting this displacement
///           when calculating star coordinates does not lead to a
///           14.6 mas change in right ascension, only a small second-
///           order distortion in the pattern of the precession-nutation
///           effect.
///
/// For these reasons, the SOFA functions do not generate the "total
/// nutations" directly, though they can of course easily be
/// generated by calling iauBi00, iauPr00 and the present function
/// and adding the results.
///
/// The MHB2000 model contains 41 instances where the same frequency
/// appears multiple times, of which 38 are duplicates and three are
/// triplicates.  To keep the present code close to the original MHB
/// algorithm, this small inefficiency has not been corrected.
///
/// @param[in] date1 (date2) TT as a 2-part Julian Date. The TT date 
///            date1+date2 is a Julian Date, apportioned in any convenient way 
///            between the two arguments. Optimally, the 'J2000 method' is best 
///            matched to the way the argument is handled internally and will 
///            deliver the optimum resolution.
/// @param[in] date2 (date1) TT as a 2-part Julian Date. See above
/// @param[out] dpsi nutation (luni-solar+planetary) longtitude component in 
///             radians
/// @param[out] deps nutation (luni-solar+planetary) obliquity component in 
///             radians
void nut00a(double date1, double date2, double &dpsi,
                               double &deps) noexcept;

/// @brief Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).
/// Naming the following points:
///       e = J2000.0 ecliptic pole,
///       p = GCRS pole,
///       E = mean ecliptic pole of date,
/// and   P = mean pole of date,
///
/// the four Fukushima-Williams angles are as follows:
///       gamb = gamma_bar = epE
///       phib = phi_bar = pE
///       psib = psi_bar = pEP
///       epsa = epsilon_A = EP
/// The matrix representing the combined effects of frame bias and
/// precession is:
///        PxB = R_1(-epsa).R_3(-psib).R_1(phib).R_3(gamb)
///
/// The matrix representing the combined effects of frame bias,
/// precession and nutation is simply:
///        NxPxB = R_1(-epsa-dE).R_3(-psib-dP).R_1(phib).R_3(gamb)
///  where dP and dE are the nutation components with respect to the
///  ecliptic of date.
/// @param[in] date1 (date2) TT as a 2-part Julian Date. The TT date 
///            date1+date2 is a Julian Date, apportioned in any convenient way 
///            between the two arguments. Optimally, The 'J2000 method' is best 
///            matched to the way the argument is handled internally and will 
///            deliver the optimum resolution.
/// @param[in] gamb F-W angle gamma_bar (radians)
/// @param[in] phib F-W angle phi_bar (radians)
/// @param[in] psib F-W angle psi_bar (radians)
/// @param[in] epsa F-W angle epsilon_A (radians)
void pfw06(double date1, double date2,
              double &gamb, double &phib, double &psib, double &epsa) noexcept;

/// @brief Form the matrix of polar motion for a given date, IAU 2000.
/// The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning that it is 
/// the final rotation when computing the pointing direction to a celestial 
/// source.
/// See iers2010, 5.4.1
/// @param[in] xp X-coordinate of the pole (radians). The arguments xp and yp 
///            are the coordinates (in radians) of the Celestial Intermediate 
///            Pole with respect to the International Terrestrial Reference 
///            System (see IERS Conventions 2003), measured along the meridians 
///            0 and 90 deg west respectively.
/// @param[in] yp Y-coordinate of the pole (radians). See above
/// @return rpom polar-motion matrix
iers2010::RotationMatrix3 pom00(double xp, double yp, double sp) noexcept;

/// @brief Compute TIO locator s'
/// The TIO locator s', positioning the Terrestrial Intermediate Origin
/// on the equator of the Celestial Intermediate Pole.
/// The TIO locator s' is obtained from polar motion observations by numerical 
/// integration, and so is in essence unpredictable. However, it is dominated 
/// by a secular drift of about 47 microarcseconds per century, which is the 
/// approximation evaluated by the present function.
/// See ier2010, 5.5.2
/// @param[in] date1  TT as a 2-part Julian Date. The TT date date1+date2 is a 
///            Julian Date, apportioned in any convenient way between the two 
///            arguments. The J2000 method is best matched to the way the 
///            argument is handled internally and will deliver the optimum 
///            resolution.
/// @param[in] date2  TT as a 2-part Julian Date (see above)
/// @return the TIO locator s' in radians
double sp00(double date1, double date2) noexcept;

/// @brief Earth rotation angle (IAU 2000 model).
/// @param[in] dj1 (dj2) UT1 as a 2-part Julian Date. The UT1 date dj1+dj2 is a 
///            Julian Date, apportioned in any convenient way between the 
///            arguments dj1 and dj2. Optimally use the 'date & time method' 
///            method.
/// @param[in] dj2 (dj1) UT1 as a 2-part Julian Date
/// @return Earth rotation angle (radians), range 0-2pi
double era00(double dj1, double dj2) noexcept;

/// @brief Compute the GCRS-to-CIRS matrix
/// Form the celestial to intermediate-frame-of-date matrix given the CIP X,Y 
/// and the CIO locator s.
/// The matrix rc2i is the first stage in the transformation from celestial to 
/// terrestrial coordinates:
/// [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
///       = RC2T * [CRS]
/// where [CRS] is a vector in the Geocentric Celestial Reference System and 
/// [TRS] is a vector in the International Terrestrial Reference System (see 
/// IERS Conventions 2003), ERA is the Earth Rotation Angle and RPOM is the 
/// polar motion matrix.
/// See iers2010, 5.4.4
/// @param[in] x Celestial Intermediate Pole, X coordinate. The Celestial 
///              Intermediate Pole coordinates are the x,y components of the 
///              unit vector in the Geocentric Celestial Reference System.
/// @param[in] y Celestial Intermediate Pole, Y coordinate (see above)
/// @param[in] s the CIO locator s. The CIO locator s (in radians) positions 
///              the Celestial Intermediate Origin on the equator of the CIP.
/// @param[out] rc2i celestial-to-intermediate matrix
void c2ixys(double x, double y, double s,
            iers2010::RotationMatrix3 &rc2i) noexcept;

/// @brief X,Y coordinates of celestial intermediate pole from series based
/// on IAU 2006 precession and IAU 2000A nutation.
/// The X,Y coordinates are those of the unit vector towards the celestial
/// intermediate pole.  They represent the combined effects of frame bias,
/// precession and nutation.
/// This routine is used for the so called 'CIO-based' transformation, see
/// iers2010, 5.5.4
/// @param[in] date1 TT as a 2-part Julian Date. The TT date date1+date2 is a
///            Julian Date, apportioned in any convenient way between the two
///            arguments. Optimally, the 'J2000' method is best matched to the 
///            way the argument is handled internally and will deliver the
///            optimum resolution.
/// @param[in] date2 (see above)
/// @param[out] x CIP X coordinate.
/// @param[out] y CIP Y coordinate
/// @note function is adopted from IAU SOFA, release 2021-05-12
void xy06(double date1, double date2, double &x, double &y) noexcept;

/// The CIO locator s, positioning the Celestial Intermediate Origin on the
/// equator of the Celestial Intermediate Pole, given the CIP's X,Y coordinates.
/// Compatible with IAU 2006/2000A precession-nutation.
/// The CIO locator s is the difference between the right ascensions of the same
/// point in two systems:  the two systems are the GCRS and the CIP,CIO, and the
/// point is the ascending node of the CIP equator.  The quantity s remains
/// below 0.1 arcsecond throughout 1900-2100.
/// The series used to compute s is in fact for s+XY/2, where X and Y are the x
/// and y components of the CIP unit vector;  this series is more compact than a
/// direct series for s would be.  This function requires X,Y to be supplied by
/// the caller, who is responsible for providing values that are consistent with
/// the supplied date.
/// The model is consistent with the "P03" precession (Capitaine et al. 2003),
/// adopted by IAU 2006 Resolution 1, 2006, and the IAU 2000A nutation (with P03
/// adjustments).
/// @param[in] date1 TT as a 2-part Julian Date. The TT date date1+date2 is a
///            Julian Date, apportioned in any convenient way between the two
///            arguments. Optimally, the 'J2000 method' is best matched to the 
///            way the argument is handled internally and will deliver the
///            optimum resolution.
/// @param[in] date2 (see above)
/// @param[in] x CIP X coordinate (see xy06)
/// @param[in] x CIP Y coordinate (see xy06)
/// @return the CIO locator s in radians
double s06(double date1, double date2, double x,
                          double y) noexcept;

/// @brief Mean obliquity of the ecliptic, IAU 2006 precession model.
/// The result is the angle between the ecliptic and mean equator of date 
/// date1+date2.
/// @param[in] date1 (date2)  TT as a 2-part Julian Date. The TT date 
///            date1+date2 is a Julian Date, apportioned in any convenient way 
///            between the two arguments. Optimally, The 'J2000 method' is best 
///            matched to the way the argument is handled internally and will 
///            deliver the optimum resolution.
/// @param[in] date2 (date1)  TT as a 2-part Julian Date. 
double obl06a(double date1, double date2) noexcept {
  // Interval between fundamental date J2000.0 and given date (JC).
  const double t = ((date1 - iers2010::DJ00) + date2) / iers2010::DJC;
  // Mean obliguity
  const double eps0 = (84381.406e0     +
          (-46.836769e0    +
          ( -0.0001831e0   +
          (  0.00200340e0  +
          ( -0.000000576e0 +
          ( -0.0000000434e0) * t) * t) * t) * t) * t) * iers2010::DAS2R;
  return eps0;
}

    /// @brief Fundamental argument, IERS Conventions (2003): mean anomaly of
    /// the
    ///        Moon.
    /// @note  - Though t is strictly TDB, it is usually more convenient to use
    /// TT,
    ///         which makes no significant difference.
    ///        - The expression used is as adopted in IERS Conventions (2003)
    ///        and is
    ///          from Simon et al. (1994).
    ///        - function is adopted from IAU SOFA, release 2021-05-12
    /// @param[in] t TDB, Julian centuries since J2000.0
    /// @return mean anomaly of the Moon in radians
    double fal03(double t) noexcept {
  double a =
      std::fmod(485868.249036 +
                    t * (1717915923.2178 +
                         t * (31.8792 + t * (0.051635 + t * (-0.00024470)))),
                iers2010::TURNAS) *
      iers2010::DAS2R;
  return a;
}

/// @brief Fundamental argument, IERS Conventions (2003): mean anomaly of the 
///        Sun.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT, 
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is 
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean anomaly of the Sun in radians
inline double falp03(double t) noexcept {
  double a =
      std::fmod(1287104.793048 +
                    t * (129596581.0481 +
                         t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))),
                iers2010::TURNAS) *
      iers2010::DAS2R;
  return a;
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of the 
///        Moon minus mean longitude of the ascending node
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT, 
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is 
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return F in radians
inline double faf03(double t) noexcept {
  double a =
      std::fmod(335779.526232 +
                    t * (1739527262.8478 +
                         t * (-12.7512 + t * (-0.001037 + t * (0.00000417)))),
                iers2010::TURNAS) *
      iers2010::DAS2R;
  return a;
}

/// @brief Fundamental argument, IERS Conventions (2003): mean elongation of 
///        the Moon from the Sun.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT, 
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is 
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return  D in radians
inline double fad03(double t) noexcept {
  double a = std::fmod(          1072260.703692 +
             t * ( 1602961601.2090 +
             t * (        - 6.3706 +
             t * (          0.006593 +
             t * (        - 0.00003169 ) ) ) ), iers2010::TURNAS ) * iers2010::DAS2R;
  return a;
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of the 
///        Moon's ascending node.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT, 
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is 
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return Omega, radians
inline double faom03(double t) noexcept {
double a = std::fmod(          450160.398036 +
             t * ( - 6962890.5431 +
             t * (         7.4722 +
             t * (         0.007702 +
             t * (       - 0.00005939 ) ) ) ), iers2010::TURNAS ) * iers2010::DAS2R;
  return a;
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of 
///        Mercury.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT, 
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is 
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean longitude of Mercury, radians
inline double fame03(double t) noexcept {
  return std::fmod(4.402608842 + 2608.7903141574 * t, iers2010::D2PI);
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of
///        Venus.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT,
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean longitude of Venus, radians
inline double fave03(double t) noexcept {
  return std::fmod(3.176146697 + 1021.3285546211 * t, iers2010::D2PI);
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of
/// Earth.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT,
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean longitude of Earth, radians
inline double fae03(double t) noexcept {
  return std::fmod(1.753470314 + 628.3075849991 * t, iers2010::D2PI);
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of
/// Mars.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT,
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean longitude of Mars, radians
inline double fama03(double t) noexcept {
  return std::fmod(6.203480913 + 334.0612426700 * t, iers2010::D2PI);
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of
/// Jupiter
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT,
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean longitude of Jupiter, radians
inline double faju03(double t) noexcept {
  return std::fmod(0.599546497 + 52.9690962641 * t, iers2010::D2PI);
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of
/// Saturn
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT,
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean longitude of Saturn, radians
inline double fasa03(double t) noexcept {
  return std::fmod(0.874016757 + 21.3299104960 * t, iers2010::D2PI);
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of
/// Uranus
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT,
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean longitude of Uranus, radians
inline double faur03(double t) noexcept {
  return std::fmod(5.481293872 + 7.4781598567 * t, iers2010::D2PI);
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of
/// Neptune
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT,
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean longitude of Neptune, radians
inline double fane03(double t) noexcept {
  return std::fmod(5.311886287 + 3.8133035638 * t, iers2010::D2PI);
}

/// @brief Fundamental argument, IERS Conventions (2003): general accumulated
/// precession in longitude.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT,
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is
///          taken from Kinoshita & Souchay (1990) and comes originally from
///          Lieske et al. (1977).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return general precession in longitude, radians
inline double fapa03(double t) noexcept {
  return (0.024381750 + 0.00000538691 * t) * t;
}

}// sofa
}// iers2010
#endif
