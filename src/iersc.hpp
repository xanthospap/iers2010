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
constexpr const double DJ00 = 2'451'545e0;

/// @brief Days per Julian century
constexpr const double DJC = 36'525e0;

/// @brief Julian Date of Modified Julian Date zero
constexpr const double DJM0 = 2'400'000.5e0;

/// @brief Reference epoch (J2000.0), Modified Julian Date
constexpr const double DJM00 = 51'544.5e0;

/// @brief Equatorial radius of the Earth [m].
constexpr const double Re = 6'378'136.6e0;

/// @brief Geocentric gravitational constant [m^3 / s ^2)]
constexpr const double GMe = 3.986004418e14;

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
