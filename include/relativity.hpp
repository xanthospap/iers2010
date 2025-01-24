/** @file
 * Relativistic effects to account for in precise geodesy and orbit 
 * determination of near-Earth orbiting satellites.
 */
#ifndef __DSO_RELATIVISTIC_CONSIDERATIONS_GEODESY_HPP__
#define __DSO_RELATIVISTIC_CONSIDERATIONS_GEODESY_HPP__

#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "iersconst.hpp"

namespace dso {

/** Relativistic correction to the acceleration of a satellite, IERS 2010.
 *
 * The computation here follows the IERS 2010 standards, see Chapter 10.3 
 * Equations of motion for an artiﬁcial Earth satellite.
 * The acceleration returned is the sum of the terms:
 * 1. The Schwarzschild terms, 
 * 2. the eﬀects of Lense-Thirring precession (frame-dragging), and
 * 3. the geodesic (de Sitter) precession
 *
 * @param[in] rsat GCRF state vector (position+velocity) of satellite in 
 *                 [m] and [m/sec] respectively.
 * @param[in] rsun GCRF state vector of Sun in [m] and [m/sec].
 * @param[in] GMe  Gravitational constant of Earth in [m^3/s^2].
 * @param[in] GMs  Gravitational constant of Sun in [m^3/s^2].
 * @param[in] J    Earth's angular momentum per unit mass in [m^2/s].
 * @param[in] c    Speed of light (vacum) in [m/sec]
 * @param[in] beta PPN (parameterized post-Newtonian) parameter β; equal to 1
 *                 for general relativity.
 * @param[in] gamma PPN (parameterized post-Newtonian) parameter γ; equal to 1
 *                 for general relativity.
 * @return         Cartesian acceleration vector in GCRF [m/sec^2]. 
 */
Eigen::Matrix<double, 3, 1> iers2010_relativistic_acceleration(
    const Eigen::Matrix<double, 6, 1> &rsat,
    const Eigen::Matrix<double, 6, 1> &rsun, double GMe = ::iers2010::GMe,
    double GMs = ::iers2010::GMs,
    double J = ::iers2010::Je, double c = ::iers2010::C, double beta = 1e0,
    double gamma = 1e0) noexcept;

} /* namespace dso */
#endif
