#ifndef __IAU_IERS10_SOFA_CPP_HPP__
#define __IAU_IERS10_SOFA_CPP_HPP__

#include "iersc.hpp"
#include <cmath>

namespace iers2010 {
  struct RotationMatrix3 {
    double data[3][3] = {{1e0,0e0,0e0}, {0e0, 1e0, 0e0}, {0e0, 0e0, 1e0}};
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
  };// RotationMatrix3
}// iers2010

namespace iers2010::sofa {

/// @brief X,Y coordinates of celestial intermediate pole from series based
/// on IAU 2006 precession and IAU 2000A nutation.
/// The X,Y coordinates are those of the unit vector towards the celestial
/// intermediate pole.  They represent the combined effects of frame bias,
/// precession and nutation.
/// This routine is used for the so called 'CIO-based' transformation, see
/// iers2010, 5.5.4
/// @param[in] date1 TT as a 2-part Julian Date. The TT date date1+date2 is a
///            Julian Date, apportioned in any convenient way between the two
///            arguments. Optimally, the first part (date1) should be the
///            reference epoch (J2000.0), Julian Date (aka 2451545.0) and date2
///            should be the remaining part of the date.
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
///            arguments. Optimally, the first part (date1) should be the
///            reference epoch (J2000.0), Julian Date (aka 2451545.0) and date2
///            should be the remaining part of the date.
/// @param[in] date2 (see above)
/// @param[in] x CIP X coordinate (see xy06)
/// @param[in] x CIP Y coordinate (see xy06)
/// @return the CIO locator s in radians
double s06(double date1, double date2, double x,
                          double y) noexcept;

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
                TURNAS) *
      DAS2R;
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
double falp03(double t) noexcept {
  double a =
      std::fmod(1287104.793048 +
                    t * (129596581.0481 +
                         t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))),
                TURNAS) *
      DAS2R;
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
double faf03(double t) noexcept {
  double a =
      std::fmod(335779.526232 +
                    t * (1739527262.8478 +
                         t * (-12.7512 + t * (-0.001037 + t * (0.00000417)))),
                TURNAS) *
      DAS2R;
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
double fad03(double t) noexcept {
  double a = std::fmod(          1072260.703692 +
             t * ( 1602961601.2090 +
             t * (        - 6.3706 +
             t * (          0.006593 +
             t * (        - 0.00003169 ) ) ) ), TURNAS ) * DAS2R;
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
double faom03(double t) noexcept {
double a = std::fmod(          450160.398036 +
             t * ( - 6962890.5431 +
             t * (         7.4722 +
             t * (         0.007702 +
             t * (       - 0.00005939 ) ) ) ), TURNAS ) * DAS2R;
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
double fame03(double t) noexcept {
  return std::fmod(4.402608842 + 2608.7903141574 * t, D2PI);
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
double fave03(double t) noexcept {
  return std::fmod(3.176146697 + 1021.3285546211 * t, D2PI);
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
double fae03(double t) noexcept {
  return std::fmod(1.753470314 + 628.3075849991 * t, D2PI);
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
double fama03(double t) noexcept {
  return std::fmod(6.203480913 + 334.0612426700 * t, D2PI);
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
double faju03(double t) noexcept {
  return std::fmod(0.599546497 + 52.9690962641 * t, D2PI);
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
double fasa03(double t) noexcept {
  return std::fmod(0.874016757 + 21.3299104960 * t, D2PI);
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
double faur03(double t) noexcept {
  return std::fmod(5.481293872 + 7.4781598567 * t, D2PI);
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
double fane03(double t) noexcept {
  return std::fmod(5.311886287 + 3.8133035638 * t, D2PI);
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
double fapa03(double t) noexcept {
  return (0.024381750 + 0.00000538691 * t) * t;
}

}// ngpt::sofa

#endif