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
#include "geodesy/geodesy.hpp"
#include "eigen3/Eigen/Eigen"

namespace dso {

class SolidEarthTide {
public:
  static constexpr const int degree = 4;

  struct PointSphericalTrigs {
    /* spherical coordinates */
    dso::SphericalCrd msph;
    /* sin (geocentric latitude) */
    double mslat;
    /* cos (geocentric latitude) */
    double mclat;
    /* sin (longitude) */
    double mslon;
    /* cos (longitude) */
    double mclon;

    template <typename T,
              typename = std::enable_if_t<CoordinateTypeTraits<T>::isCartesian>>
    explicit PointSphericalTrigs(const /*dso::CartesianCrd*/ T &rsta) noexcept
        : msph(dso::cartesian2spherical(rsta)) {
      mslat = std::sin(msph.lat());
      mclat = std::cos(msph.lat());
      mslon = std::sin(msph.lon());
      mclon = std::cos(msph.lon());
    }
  }; /* struct PointsSphericalTrigs */

private:
  /* gravitational constant of Sun/Earth */
  double mGMSun;
  /* gravitational constant of Moon/Earth */
  double mGMMoon;
  /* Mass of Moon / Mass of Earth */
  double mmeratio{dso::MoonEarthMassRatio};
  /* Mass of Moon / Mass of Earth */
  double mseratio{332946.0482e0};
  /* Stokes coefficients of (n,m) = (4,4) */
  StokesCoeffs mcs;

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
  int potential_step1(const Eigen::Matrix<double, 3, 1> &rMoon,
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
  int potential_step2(const MjdEpoch &mjdtt, const MjdEpoch &mjdut1,
                             const double *const delaunay_args, double &dC20,
                             double &dC21, double &dS21, double &dC22,
                             double &dS22) const noexcept;

  /** @brief Compute Step-1 displacement due to solid earth tides, according to
   *   IERS 2010.
   *
   * Step-1 corrections, include:
   *  + Displacement due to degree n = 2,3 tide (in-phase).
   *    The formulation here follows the IERS 2010 standrads, Section 7.1.1.1
   *    Conventional model for solid Earth tides.
   *    This function treats the Step-1 in-phase corrections for degrees 
   *    n=2,3 following Equations (5) and (6) respectively.
   *
   *  + Displacement due to degree n=2 tide (out-of-phase) and the contribution
   *    from latitude dependence.
   *    The formulation here follows the IERS 2010 standrads, Section 7.1.1.1
   *    Conventional model for solid Earth tides.
   *    This function treats the Step-1 out-of-phase corrections for degree 
   *    n=2, following Equations (10) and (11) for the diurnal and 
   *    semi-diurnal tides respectively.
   * 
   *  + Additionally, the function add the contribution from latitude 
   *    dependence (Equations (8) and (9)) for the diurnal and semi-diurnal 
   *    band.
   *
   * @param[in] rsta  ECEF, cartesian position vector of site/point on Earth
   * @param[in] rsun  ECEF, cartesian position vector of Sun
   * @param[in] rmoon ECEF, cartesian position vector of Moon
   * @param[in] tsta  An instance of type SolidEarthTide::PointSphericalTrigs
   *            holding spherical coordinates and trigonometric numbers of
   *            the site
   * @param[in] tsun  An instance of type SolidEarthTide::PointSphericalTrigs
   *            holding spherical coordinates and trigonometric numbers of
   *            the Sun
   * @param[in] tmoon An instance of type SolidEarthTide::PointSphericalTrigs
   *            holding spherical coordinates and trigonometric numbers of
   *            the Moon
   * @return Displacement vector in local tangent coordinates, where the
   *         Up is in the radial (not vertical) direction and East and North
   *         are perpendicular to that direction, i.e.
   *         Δr = [East, North, Up/Radial]
   */
  Eigen::Matrix<double, 3, 1> step1_displacement(
      const Eigen::Matrix<double, 3, 1> &rsta,
      const Eigen::Matrix<double, 3, 1> &rMoon,
      const Eigen::Matrix<double, 3, 1> &rSun,
      const dso::SolidEarthTide::PointSphericalTrigs &tSta,
      const dso::SolidEarthTide::PointSphericalTrigs &tMoon,
      const dso::SolidEarthTide::PointSphericalTrigs &tSun) noexcept;

  /** This function computes corrections to Step-1 Solid Earth tide -induced 
   * displacement at a given point, due to frequency dependence of Love 
   * and Shida numbers.
   * Both diurnal (for degree n=2) and long-period tides (n=2) are considered 
   * here.
   *
   * The model follows IERS 2010 Conventions, Section 7.1.1.1 Conventional 
   * model for solid Earth tides. Equations involved are (mainly) Eq. 12 and 
   * Eq. 13, while the consodered frequencies are taken from Tables 7.3a and 
   * 7.3b.
   *
   * @param[in] mjdtt  Time/epoch of computation in [TT]
   * @param[in] mjdut1 Time/epoch of computation in [UT1]
   * @param[in] rsta ECEF, cartesian coordinates of site/point on Earth [m]
   * @param[in] tSta An instance of type PointSphericalTrigs, holding trig 
   *            numbers for rsta
   * @param[in] delaunay_args Fundamental/Delaunay arguments at the time of 
   *            computation, i.e. [l, lp, f, d, Ω] in [rad]
   * @return Displacement vector in local tangent coordinates, where the
   *         Up is in the radial (not vertical) direction and East and North
   *         are perpendicular to that direction, i.e.
   *         Δr = [East, North, Up/Radial]
   */
  Eigen::Matrix<double, 3, 1>
  step2_displacement(const MjdEpoch &mjdtt, const MjdEpoch &mjdut1,
                     const Eigen::Matrix<double, 3, 1> &rsta,
                     const dso::SolidEarthTide::PointSphericalTrigs &tsta,
                     const double *const fargs) noexcept;

public:
  /// @brief Constructor
  /// @param GMearth Standard gravitational parameter μ=GM for the Earth
  /// [m^2/s^2]
  /// @param Rearth Equatorial radius of Earth [m]
  /// @param GMmoon Standard gravitational parameter μ=GM for the Moon [m^2/s^2]
  /// @param GMsun Standard gravitational parameter μ=GM for the Moon [m^2/s^2]
  SolidEarthTide(double GMearth, double Rearth, double GMmoon,
                 double GMsun) noexcept
      : mGMSun(GMsun), mGMMoon(GMmoon), mcs(degree, degree, GMearth, Rearth){};

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
  
  /** Compute site displacement due to Solid Earth tide according to IERS 2010.
   *
   * This function follows the IERS 2010 model, using a two-step approach.
   * Both Step-1 (due to Sun and Moon) and Step-2 (frequency-dependent) 
   * corrections are considered.
   *
   * @param[in] mjdtt  Time/epoch of computation in [TT]
   * @param[in] mjdut1 Time/epoch of computation in [UT1]
   * @param[in] rMoon ECEF coordinates of moon [m]
   * @param[in] rSun  ECEF coordinates of Sun [m]
   * @param[in] delaunay_args Fundamental/Delaunay arguments at the time of 
   *            computation, i.e. [l, lp, f, d, Ω] in [rad]
   * @return Always 0
   */
  Eigen::Matrix<double, 3, 1>
  displacement(const MjdEpoch &mjdtt, const MjdEpoch &mjdut1,
               const Eigen::Matrix<double, 3, 1> &rsta,
               const Eigen::Matrix<double, 3, 1> &rMoon,
               const Eigen::Matrix<double, 3, 1> &rSun,
               const double *const delaunay_args) noexcept;
}; /* SolidEarthTide */
} /* namespace dso */

#endif
