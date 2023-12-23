#include "solid_earth_tide.hpp"
#include "geodesy/geodesy.hpp"
#include <cmath>
#include <cstdio>
#include <geodesy/crd_transformations.hpp>

namespace {
/** Nominal Love & Shida numbers */
constexpr const double love_h0 = 0.6078e0;
constexpr const double love_h2 = -0.0006e0;
constexpr const double love_h3 = 0.292e0;
constexpr const double shida_l0 = 0.0847e0;
constexpr const double shida_l2 = 0.0002e0;
constexpr const double shida_l3 = 0.015e0;

/** @brief Compute displacement due to degree n = 2,3 tide (in-phase).
 *
 * The formulation here follows the IERS 2010 standrads, Section 7.1.1.1
 * Conventional model for solid Earth tides.
 * This function treats the Step-1 in-phase corrections for degrees n=2,3,
 * following Equations (5) and (6) respectively.
 *
 * @param[in] Re Equatorial radius of the Earth [m]
 * @param[in] GM Gravitational constant of Earth
 * @param[in] rtb ECEF, cartesian coordinates of third body (sun or
 *            Moon) [rad] and [m]
 * @param[in] GMtb Gravitational constant of Third body
 * @param[in] rsta ECEF, spherical coordinates of site/point on Earth
 * @return The displacement vector Δr in ECEF Cartesian coordinates, i.e. 
 *         Δr = [δX, δY, δZ]
 */
Eigen::Matrix<double, 3, 1> step1_in_phase(
    double Re, double GM, const Eigen::Matrix<double, 3, 1> &rtb, double GMtb,
    const Eigen::Matrix<double, 3, 1> &rsta,
    const dso::SolidEarthTide::PointSphericalTrigs &statrigs) noexcept {

  /* degree 2 Love number, taking latitude into consideration, Eq. (2) */
  const double __sf = statrigs.mslat;
  const double __3sfm1d2 = (3e0 * __sf * __sf - 1e0) / 2e0;
  const double love2 = love_h0 + love_h2 * __3sfm1d2;
  /* degree 2 Shida number, taking latitude into consideration, Eq. (2) */
  const double shida2 = shida_l0 + shida_l2 * __3sfm1d2;

  /* displacement vector (to be returned) */
  Eigen::Matrix<double, 3, 1> dr = Eigen::Matrix<double, 3, 1>::Zero();

  /* unit vector from the geocenter to third body */
  const auto urtb = rtb.normalized();
  /* unit vector from the geocenter to site */
  const auto ursta = rsta.normalized();
  /* Eq. (5), displacement due to degree n=2 */
  const auto __Rjr = urtb.dot(ursta);
  dr = love2 * ursta * ((3e0 * __Rjr * __Rjr - 1e0) / 2e0) +
       3e0 * shida2 * __Rjr * (urtb - __Rjr * ursta);
  /* Eq. (6), displacement due to degree n=3 */
  const double scale3 = Re / rtb.norm();
  dr += scale3 *
        (love_h3 * ursta *
             ((5e0 / 2e0) * __Rjr * __Rjr - (3e0 / 2e0) * __Rjr) +
         shida_l3 * ((15e0 / 2e0) * __Rjr * __Rjr -
                     (3e0 / 2e0)) * (urtb - __Rjr * ursta));

  /* scale and return */
  const double scale = (GMtb / GM) * std::pow(Re / rtb.norm(), 3) * Re;
  return scale * dr;
}

/** @brief Compute displacement due to degree n=2 tide (out-of-phase) and the
 *   constribution from latitude dependence.
 *
 * The formulation here follows the IERS 2010 standrads, Section 7.1.1.1
 * Conventional model for solid Earth tides.
 * This function treats the Step-1 out-of-phase corrections for degree n=2,
 * following Equations (10) and (11) for the diurnal and semi-diurnal
 * tides respectively.
 *
 * Additionally, the function add the contribution from latitude dependence
 * (Equations (8) and (9)) for the diurnal and semi-diurnal band.
 *
 * @param[in] Re Equatorial radius of the Earth [m]
 * @param[in] GM Gravitational constant of Earth
 * @param[in] rtb ECEF, cartesian position vector of third body (sun or
 *            Moon) [m]
 * @param[in] GMtb Gravitational constant of Third body
 * @param[in] rsta ECEF, cartesian position vector of site/point on Earth
 * @return Displacement vector in local tangent coordinates, where the 
 *         Up is in the radial (not vertical) direction and East and North 
 *         are perpendicular to that direction, i.e.
 *         Δr = [East, North, Up/Radial]
 */
Eigen::Matrix<double, 3, 1> step1_outof_phase(
    double Re, double GM, const Eigen::Matrix<double, 3, 1> &rtb, double GMtb,
    [[maybe_unused]]const Eigen::Matrix<double, 3, 1> &rsta,
    const dso::SolidEarthTide::PointSphericalTrigs &tbtrigs,
    const dso::SolidEarthTide::PointSphericalTrigs &statrigs) noexcept {

  /* var trigonometric numbers */
  const double __s2Phi = 2e0 * tbtrigs.mslat * tbtrigs.mclat;
  const double __cPhi = tbtrigs.mclat;
  const double __s2phi = 2e0 * statrigs.mslat * statrigs.mclat;
  const double __c2phi =
      statrigs.mclat * statrigs.mclat - statrigs.mslat * statrigs.mslat;
  const double __sphi = statrigs.mslat;
  const double __cphi = statrigs.mclat;
  const double __sPhi = tbtrigs.mslat;
  const double __slmL = std::sin(statrigs.msph.lon() - tbtrigs.msph.lon());
  const double __clmL = std::cos(statrigs.msph.lon() - tbtrigs.msph.lon());
  const double __s2lmL = 2e0 * __slmL * __clmL;
  const double __c2lmL = __clmL * __clmL - __slmL * __slmL;

  /* displacement vector (to be returned) */
  Eigen::Matrix<double, 3, 1> dr = Eigen::Matrix<double, 3, 1>::Zero();

  /* Eq. (10a) for radial and (10b) traverse directions; diurnal
   * tides (n=2)
   */
  {
    const double love_hI10 = -0.0025e0;
    const double shida_lI10 = -0.0007e0;
    /* east */
    dr(0) = -(3e0 / 2e0) * shida_lI10 * __s2Phi * __sphi * __clmL;
    /* north */
    dr(1) = -(3e0 / 2e0) * shida_lI10 * __s2Phi * __c2phi * __slmL;
    /* radial */
    dr(2) = -(3e0 / 4e0) * love_hI10 * __s2Phi * __s2phi * __slmL;
  }

  /* Eq. (11a) for radial and (11b) traverse directions; semi-diurnal
   * tides (n=2)
   */
  {
    const double love_hI11 = -0.0022e0;
    const double shida_lI11 = -0.0007e0;
    const double __cPhi2 = __cPhi * __cPhi;
    const double __cphi2 = __cphi * __cphi;
    /* east */
    dr(0) += -(3e0 / 2e0) * shida_lI11 * __cPhi2 * __cphi * __c2lmL;
    /* north */
    dr(1) += (3e0 / 4e0) * shida_lI11 * __cPhi2 * __s2phi * __s2lmL;
    /* radial */
    dr(2) += -(3e0 / 4e0) * love_hI11 * __cPhi2 * __cphi2 * __s2lmL;
  }

  /* contribution from the diurnal band, with l(1) = 0.0012, i.e.
   * from latitude dependence. Eq. (8)
   */
  {
    const double shida_l18 = 0.0012e0;
    const double P21 = -3e0 * __sPhi * __cPhi;
    /* east */
    dr(0) += shida_l18 * __sphi * P21 * __c2phi * __slmL;
    /* north */
    dr(1) -= shida_l18 * __sphi * P21 * __sphi * __clmL;
  }

  /* contribution from the semi-diurnal band, with l(1) = 0.0024, i.e.
   * from latitude dependence. Eq. (9)
   */
  {
    const double shida_l19 = 0.0024e0;
    const double P22 = 3e0 * __cPhi * __cPhi;
    /* east */
    dr(0) += -.5e0 * shida_l19 * __sphi * __cphi * P22 * __sphi * __s2lmL;
    /* north */
    dr(1) += -.5e0 * shida_l19 * __sphi * __cphi * P22 * __c2lmL;
  }

  /* scale and return */
  const double scale = (GMtb / GM) * std::pow(Re / rtb.norm(), 3) * Re;
  return scale * dr;
}
} /* unnamed namespace */

Eigen::Matrix<double, 3, 1> dso::SolidEarthTide::step1_displacement(
    const Eigen::Matrix<double, 3, 1> &rsta,
    const Eigen::Matrix<double, 3, 1> &rMoon,
    const Eigen::Matrix<double, 3, 1> &rSun,
    const dso::SolidEarthTide::PointSphericalTrigs &tSta,
    const dso::SolidEarthTide::PointSphericalTrigs &tMoon,
    const dso::SolidEarthTide::PointSphericalTrigs &tSun) noexcept {

  /* step-1 corrections (time domain) for Sun (Cartesian, ECEF) */
  Eigen::Matrix<double, 3, 1> drxyz =
      step1_in_phase(mcs.Re(), mcs.GM(), rSun, mGMSun, rsta, tSta);
  drxyz += step1_in_phase(mcs.Re(), mcs.GM(), rMoon, mGMMoon, rsta, tSta);

  /* get rotation matrix to transform between Cartesian and topocentric, i.e.
   * [x,y,z] = R *[enu]
   */
  const auto R = dso::geodetic2lvlh(tSta.msph.lat(), tSta.msph.lon());
  
  /* transform displacement vector from Cartesian to topocentric (enu) */
  Eigen::Matrix<double, 3, 1> dr = R.transpose() * drxyz;

  /* step1 : out-of-phase (+lat.dependence) */
  dr += step1_outof_phase(mcs.Re(), mcs.GM(), rSun, mGMSun, rsta, tSun, tSta);
  dr +=
      step1_outof_phase(mcs.Re(), mcs.GM(), rMoon, mGMMoon, rsta, tMoon, tSta);
  printf("\tStep1 :%+.6f %+.6f %+.6f (Sun+Moon)\n", dr(0), dr(1), dr(2));

  /* all done */
  return dr;
}
