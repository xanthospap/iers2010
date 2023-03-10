#ifndef __IERS2010_GEN_CONSTANTS_HPP__
#define __IERS2010_GEN_CONSTANTS_HPP__

/// @file This file includes numeric constants defined in IERS2010 standards
/// @see IERS Technical Note N.36, Chapter 1, Table 1.1

#include <cmath>

namespace iers2010 {

/// @brief Pi constant
constexpr const double DPI = 3.141592653589793238462643e0;

/// @brief 2Pi constant
constexpr const double D2PI = 6.283185307179586476925287e0;

/// @brief Arcseconds in a full circle
constexpr const double TURNAS = 1'296'000e0;

/// @brief Arcseconds to radians
constexpr const double DAS2R = 4.848136811095359935899141e-6;

/// @brief Radians to arcseconds
constexpr double DR2AS = 206'264.8062470963551564734e0;

/// @brief Milliarcseconds to radians
constexpr const double DMAS2R = DAS2R / 1e3;

/// @brief Reference epoch (J2000.0), Julian Date
/// constexpr const double DJ00 = 2'451'545e0; dso::j2000_jd

/// @brief Days per Julian century
/// constexpr const double DJC = 36'525e0; dso::days_in_julian_cent

/// @brief Julian Date of Modified Julian Date zero
constexpr const double DJM0 = 2'400'000.5e0;

/// @brief Reference epoch (J2000.0), Modified Julian Date
constexpr const double DJM00 = 51'544.5e0;

/// @brief Equatorial radius of the Earth [m].
/// @see Table 1.1: IERS numerical standards, IERS 2010
constexpr const double Re = 6'378'136.6e0;

/// @brief Geocentric gravitational constant [m^3 / s ^2)]
/// @see Table 1.1: IERS numerical standards, IERS 2010
constexpr const double GMe = 3.986004418e14;

/// @brief Astronomical unit [m]
constexpr const double AU = 1.49597870700e11;

/// @brief Nominal mean Earth’s angular velocity [rads/sec]
/// @note Table 1.2: Parameters of the Geodetic Reference System GRS80, 
///       IERS2010
/// constexpr const double OmegaEarth = 7.292115e-5;
/// see https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html
constexpr const double OmegaEarth = 7.2921151467064e-5;


/// @brief Solar radius [m]
/// @see IAU, Resolution B3 on recommended nominal conversion constants for 
///      selected solar and planetary properties, August 13, 2015
constexpr const double Rs = 6.957e8;

/// @brief 1 - d(TT)/d(TCG)
/// @see Table 1.1: IERS numerical standards, IERS 2010
/// @note This quantity is U_geo / c^2, where U_geo is the potential on the
///      Earth's geoid.
constexpr const double LG = 6.9692901341e-10;

/// @brief Speed of light [m/sec]
/// @see Table 1.1: IERS numerical standards, IERS 2010
constexpr const double C = 299792458e0;

/// @brief Heliocentric gravitational constant in [m^2/s^2]
/// @see Table 1.1: IERS numerical standards, IERS 2010
constexpr const double GMSun = 1.32712442099e20;

/// @brief Normalize angle into the range 0 <= a < 2pi.
/// @param[in] angle Angle in radians
/// @return Angle in radians in range [0,2pi)
inline double nang_02pi(double angle) noexcept {
  double w = std::fmod(angle, iers2010::D2PI);
  if (w < 0e0)
    w += iers2010::D2PI;
  return w;
}

} // iers2010

#endif
