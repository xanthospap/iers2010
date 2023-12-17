/** @file
 * Solid Earth Tide Model & computation for Geodesy.
 * Solid Earth Tide has an effect both on the geopotential (expressible via 
 * spherical harmonics coefficients) and on points/sites on the Earth's 
 * surface, i.e. causing deformation.
 * Both effects are treated here in speration.
 */

#ifndef __DSO_SOLID_EARTH_TIDE_HPP__
#define __DSO_SOLID_EARTH_TIDE_HPP__

#include "stokes_coefficients.hpp"
#include "datetime/calendar.hpp"

namespace dso {

class SolidEarthTide {
  static constexpr const int degree = 4;

private:
  StokesCoeffs m_cs;

  /* @brief Compute the Step-1 effect of Solid Earth Tides on geopotential
   *  coefficient corrections ΔC_{nm} and ΔS_{nm}.
   *
   * Compute the Step-1 effect of Solid Earth tides, as described in IERS 2010,
   * Sec. 6.2.1 (Eq. 6.6 and Eq. 6.7). That is, with frequency-independent
   * values knm (i.e. Love numbers, see IERS 2010, Table 6.3), changes induced
   * by the (nm) part of the TGP in the normalized geopotential coeﬃcients are
   * computed as corrections (to the normalized coefficients), ie. ΔC{nm} and
   * ΔS{nm}. Affects the ΔC_nm and ΔS_nm (correction) coefficients, for:
   *
   *        (nm) = (20), (3,0), (4,0)
   *               (21), (3,1), (4,1)
   *               (22), (3,2), (4,2)
   *                     (3,3)
   *               -----|------|------
   *                6.6   6.6    6.7      [IERS2010 Equation nr]
   *
   * The 'third bodies' considered here are the Sun and Moon.
   *
   * @param[in] rMoon ECEF coordinates of moon [m]
   * @param[in] rSun  ECEF coordinates of Sun [m]
   * @param[out] dC   Normalized corrections to C coefficients, in the order:
   *             dC = C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,0,0
   * @param[out] dS   Normalized corrections to S coefficients, in the order:
   *             dS = 0,S21,S22,0,S31,S32,S33,0,S41,S42,0,0
   */
  int solid_earth_tide_step1(const Eigen::Matrix<double, 3, 1> &rMoon,
                             const Eigen::Matrix<double, 3, 1> &rSun,
                             std::array<double, 12> &dC,
                             std::array<double, 12> &dS) noexcept;

  /* @brief Compute the Step-2 effect of Solid Earth Tides on geopotential
   *  coefficient corrections ΔC_{2m} and ΔS_{2m}.
   *
   * Compute the Step-2 effect/corrections of Solid Earth tides, as described 
   * in IERS 2010, Sec. 6.2.1 (Eq. 6.8a through Eq. 6.7e). 
   * These corrections affect only the degree n=2 geopotential coefficients 
   * (i.e. C20, C21, C21, S21 and S22) and are used to apply small, frequency 
   * dependent corrections to the Step-1 geopotential coeffs (ΔC and ΔS). 
   * These are computed as the sums of contributions from a number of tidal 
   * constituents belonging to the respective bands.
   *
   * @param[in] mjdtt  Time/epoch of computation in [TT]
   * @param[in] mjdut1 Time/epoch of computation in [UT1]
   * @param[in] delaunay_args Fundamental/Delaunay arguments at the time of 
   *            computation, i.e. [l, lp, f, d, Ω] in [rad]
   * @param[out] dC20, dC21, dS21, dC22, dS22 Step-2 corrections for the 
   *            respective geopotential coefficients.
   */
  int solid_earth_tide_step2(const MjdEpoch &mjdtt, const MjdEpoch &mjdut1,
                             const double *const delaunay_args, double &dC20,
                             double &dC21, double &dS21, double &dC22,
                             double &dS22) const noexcept;

public:
  /// @brief Constructor
  /// @param GMearth Standard gravitational parameter μ=GM for the Earth
  /// [m^2/s^2]
  /// @param Rearth Equatorial radius of Earth [m]
  /// @param GMmoon Standard gravitational parameter μ=GM for the Moon [m^2/s^2]
  /// @param GMsun Standard gravitational parameter μ=GM for the Moon [m^2/s^2]
  SolidEarthTide(double GMearth, double Rearth, double GMmoon,
                 double GMsun) noexcept {};

  int operator()();
}; /* SolidEarthTide */
} /* namespace dso */

#endif
