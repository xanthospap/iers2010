/** @file 
 * This file includes numeric constants defined in IERS2010 standards
 *  @see IERS Technical Note N.36, Chapter 1, Table 1.1
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

/** @brief Heliocentric gravitational constant in [m^2/s^2]
 *  @see Table 1.1: IERS numerical standards, IERS 2010
 */
constexpr const double GMSun = 1.32712442099e20;

/* @brief Mean equatorial gravity in [m/sec^2], see
 * Table 1.1: IERS numerical standards, IERS 2010
 */
constexpr const double ge = 9.7803278e0;

} /* iers2010 */

#endif
