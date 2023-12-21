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
#include "eigen3/Eigen/Eigen"

namespace dso {

class SolidEarthTide {
  static constexpr const int degree = 4;

private:
  /* gravitational constant of Sun */
  double mGMSun;
  /* gravitational constant of Moon */
  double mGMMoon;
  /* Stokes coefficients of (n,m) = (4,4) */
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
   *        (nm) = (2,0), (3,0), (4,0)
   *               (2,1), (3,1), (4,1)
   *               (2,2), (3,2), (4,2)
   *                      (3,3)
   *               ------|------|------
   *                6.6    6.6     6.7   [IERS2010 Equation nr]
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

  /** Compute corrections to geopotemtial coefficients (ΔC, ΔS) due to Solid 
   * Earth tide according to IERS 2010.
   *
   * This function follows the IERS 2010 model, using a two-step approach.
   * Both Step-1 (due to Sun and Moon) and Step-2 (frequency-dependent) 
   * corrections are considered.
   * The degree/order corrections actually computed span n=2,3,4; for n=4, 
   * max(m)=2 (i.e. we compute (4,0), (4,1) and (4,2)). However, the instance 
   * will hold a Stokes coefficient matrix of (n,m)=(4,4), placing zero values
   * where needed.
   * After the computation, the instance's m_cnm member will be filled with 
   * the computed geopotential coefficient corrections (computed here). Hence, 
   * the m_cnm instance could be added to the values of a geopotential model 
   * to account for the Solid Earth Tide effect.
   *
   * @param[in] mjdtt  Time/epoch of computation in [TT]
   * @param[in] mjdut1 Time/epoch of computation in [UT1]
   * @param[in] rMoon ECEF coordinates of moon [m]
   * @param[in] rSun  ECEF coordinates of Sun [m]
   * @param[in] delaunay_args Fundamental/Delaunay arguments at the time of 
   *            computation, i.e. [l, lp, f, d, Ω] in [rad]
   * @return Always 0
   */
  int stokes_coeffs(const MjdEpoch &mjdtt, const MjdEpoch &mjdut1,
                    const Eigen::Matrix<double, 3, 1> &rMoon,
                    const Eigen::Matrix<double, 3, 1> &rSun,
                    const double *const delaunay_args) noexcept;
}; /* SolidEarthTide */
} /* namespace dso */

#endif
