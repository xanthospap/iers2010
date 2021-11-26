#ifndef __IERS2010_GEN_CONSTANTS_HPP__
#define __IERS2010_GEN_CONSTANTS_HPP__

namespace ier2010 {
/// @brief Arcseconds in a full circle
constexpr double TURNAS  =1296000e0;

/// @brief Arcseconds to radians
constexpr double DAS2R = 4.848136811095359935899141e-6;

/// @brief Reference epoch (J2000.0), Julian Date
constexpr double DJ00 = 2451545e0;

/// @brief Days per Julian century
constexpr double DJC = 36525e0;
}// iers2010

#endif