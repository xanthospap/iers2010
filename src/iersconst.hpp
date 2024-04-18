/** @file 
 * This file includes numeric constants defined in IERS2010 standards
 * @see IERS Technical Note N.36, Chapter 1, Table 1.1
 *
 * References:
 * [1] Pitjeva, Elena V. and Erland Myles Standish. “Proposals for the masses 
 * of the three largest asteroids, the Moon-Earth mass ratio and the 
 * Astronomical Unit.” Celestial Mechanics and Dynamical Astronomy 103 
 * (2009): 365-372.
 */

#ifndef __IERS2010_GEN_CONSTANTS_HPP__
#define __IERS2010_GEN_CONSTANTS_HPP__

namespace iers2010 {

/** @brief Equatorial radius of the Earth [m].
 *  @see Table 1.1: IERS numerical standards, IERS 2010
 */
constexpr const double Re = 6'378'136.6e0;

/** @brief Geocentric gravitational constant [m^3 / s^2]
 *  @see Table 1.1: IERS numerical standards, IERS 2010
 */
constexpr const double GMe = 3.986004418e14;

/* @brief Constant of gravitation in [m^3/kg/sec-2]
 * @see Table 1.1: IERS numerical standards, IERS 2010
 */
constexpr const double G = 6.67428e-11;

/** @brief Astronomical unit [m] */
constexpr const double AU = 1.49597870700e11;

/** @brief Nominal mean Earth’s angular velocity [rads/sec]
 * @note Table 1.2: Parameters of the Geodetic Reference System GRS80, 
 *       IERS2010
 * constexpr const double OmegaEarth = 7.292115e-5;
 * see https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html
 */
constexpr const double OmegaEarth = 7.2921151467064e-5;

/** @brief Solar radius [m]
 * @see IAU, Resolution B3 on recommended nominal conversion constants for 
 *      selected solar and planetary properties, August 13, 2015
 */
constexpr const double Rs = 6.957e8;

/** @brief 1 - d(TT)/d(TCG)
 * @see Table 1.1: IERS numerical standards, IERS 2010
 * @note This quantity is U_geo / c^2, where U_geo is the potential on the
 *      Earth's geoid.
 */
constexpr const double LG = 6.9692901341e-10;

/** @brief Speed of light [m/sec]
 *  @see Table 1.1: IERS numerical standards, IERS 2010
 */
constexpr const double C = 299792458e0;

/** @brief Heliocentric gravitational constant in [m^3/s^2]
 *  @see Table 1.1: IERS numerical standards, IERS 2010
 */
constexpr const double GMs = 1.32712442099e20;

/* @brief Mean equatorial gravity in [m/sec^2], see
 * Table 1.1: IERS numerical standards, IERS 2010
 */
constexpr const double ge = 9.7803278e0;

/** @brief Earth's angular momentum per unit mass in [m^2/s]. See IERS 2010, 
 * Section 10.3 Equations of motion for an artiﬁcial Earth satellite.
 */
constexpr const double Je = 9.8e8;

} /* iers2010 */

namespace dso {
/** Constants not included in IERS 2010 Standards. These are defined within 
 * the dso namespace
 */

/** Moon-Earth mass ratio according to [1] */
constexpr const double MoonEarthMassRatio = 0.0123000371e0;
/** Earth–Moon mass ratio according to [1] */
constexpr const double EarthMoonMassRatio = 81.300568e0;
/** Sun-Earth Mass Ratio according to [1] */
constexpr const double SunEarthMassRatio = 332946.0487e0;
/** Mass of the Sun according to [1] in [kg] */
constexpr const double SunMass =  1.9884e30;
/** Mass of the Earth according to [1] in [kg] */
constexpr const double SunEarth =  5.9722e24;
} /* namespace dso */

#endif
