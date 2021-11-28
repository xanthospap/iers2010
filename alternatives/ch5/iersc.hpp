#ifndef __IERS2010_GEN_CONSTANTS_HPP__
#define __IERS2010_GEN_CONSTANTS_HPP__

#include <cmath>
namespace iers2010 {
/// @brief Pi constant
constexpr double DPI = 3.141592653589793238462643e0;

/// @brief 2Pi constant
constexpr double D2PI = 6.283185307179586476925287e0;

/// @brief Arcseconds in a full circle
constexpr double TURNAS  =1296000e0;

/// @brief Arcseconds to radians
constexpr double DAS2R = 4.848136811095359935899141e-6;

// Milliarcseconds to radians
constexpr double DMAS2R = DAS2R / 1e3;

/// @brief Reference epoch (J2000.0), Julian Date
constexpr double DJ00 = 2451545e0;

/// @brief Days per Julian century
constexpr double DJC = 36525e0;

// Julian Date of Modified Julian Date zero
constexpr double DJM0 = 2400000.5e0;

// Reference epoch (J2000.0), Modified Julian Date
constexpr double DJM00 = 51544.5e0;

/// @brief Normalize angle into the range 0 <= a < 2pi.
/// @param[in] angle Angle in radians
/// @return Angle in radians in range [0,2pi)
inline double nang_02pi(double angle) noexcept {
    double w = std::fmod(angle, iers2010::D2PI);
    if (w<0e0) w += iers2010::D2PI;
    return w;
}
}// iers2010

#endif
