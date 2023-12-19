#include "solid_earth_tide.hpp"
#include <cmath>

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
 * @param[in] rtb ECEF, cartesian position vector of third body (sun or
 *            Moon) [m]
 * @param[in] GMtb Gravitational constant of Third body
 * @param[in] rsta ECEF, cartesian position vector of site/point on Earth
 *
 * TODO return ?
 */
Eigen::Matrix<double, 3, 1>
step1_in_phase(double Re, double GM, const Eigen::Matrix<double, 3, 1> &rtb,
               double GMtb, const Eigen::Matrix<double, 3, 1> &sta) noexcept {

  /* get geocentric latitude and longitude of site */
  // TODO phi, lambda

  /* degree 2 Love number, taking latitude into consideration, Eq. (2) */
  const double slat = std::sin(lat);
  const double __3sfm1d2 = (3e0 * slat * slat - 1e0) / 2e0;
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
       3e0 * shida_l2 * __Rjr * (urtb - __Rjr * ursta);
  /* Eq. (6), displacement due to degree n=3 */
  const double scale3 = Re / rtb.norm();
  dr += scale3 *
        (love_h3 * ursta *
             ((5e0 / 2e0) * __Rjr * __Rjr * __Rjr - (3e0 / 2e0) * __Rjr) +
         shida_l3 * ((15e0 / 2e0) * __Rjr * __Rjr -
                     (3e0 / 2e0) * (urtb - __Rjr * ursta)));

  /* scale and return */
  const double scale = (GM / GMtb) * std::pow(Re / rtb.norm(), 3) * Re;
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
 *
 * TODO return ?
 */
Eigen::Matrix<double, 3, 1>
step1_outof_phase(double Re, double GM, const Eigen::Matrix<double, 3, 1> &rtb,
                  double GMtb,
                  const Eigen::Matrix<double, 3, 1> &sta) noexcept {
  /* get geocentric latitude and longitude of third body */
  // TODO Phi, Lambda
  /* get geocentric latitude and longitude of site */
  // TODO phi lambda
  // unit vector in east direction: ue
  // unit vector perpendicular to ˆr in the northward direction. un

  const double __s2Phi = std::sin(2e0 * Phi);
  const double __s2phi = std::sin(2e0 * phi);
  const double __sphi = std::sin(phi);
  const double __sPhi = std::sin(Phi);
  const double __cPhi = std::cos(Phi);
  const double __c2phi = std::cos(2e0 * phi);
  const double __slmL = std::sin(lambda - Lambda);
  const double __clmL = std::cos(lambda - Lambda);

  /* displacement vector (to be returned) */
  Eigen::Matrix<double, 3, 1> dt = Eigen::Matrix<double, 3, 1>::Zero();
  double radial = 0e0;

  /* Eq. (10a) for radial and (10b) traverse directions; diurnal
   * tides (n=2)
   */
  {
    const double love_hI10 = -0.0025e0 const double shida_lI10 =
        -0.0007e0 radial =
            -(3e0 / 4e0) * love_hI10 * __s2Phi * __s2phi * __slmL;
    dt = -(3e0 / 2e0) * shida_lI10 * __s2Phi *
         (__c2phi * __slmL * un + __sphi * __clmL * ue);
  }

  /* Eq. (11a) for radial and (11b) traverse directions; semi-diurnal
   * tides (n=2)
   */
  {
    const double love_hI11 = -0.0022e0 const double shida_lI11 =
        -0.0007e0 double radial =
            -(3e0 / 4e0) * love_hI10 * __s2Phi * __s2phi * __slmL;
    dt += -(3e0 / 2e0) * shida_lI10 * __s2Phi *
          (__c2phi * __slmL * un + __sphi * __clmL * ue);
  }

  /* contribution from the diurnal band, with l(1) = 0.0012, i.e.
   * from latitude dependence. Eq. (8)
   */
  {
    const double shida_l18 = 0.0012e0;
    const double P21 = -3e0 * __sPhi * __cPhi;
    dt += -shida_l18 * __sphi * P21 *
          (__sphi * __clmL * un - __c2phi * __slmL * ue);
  }

  /* contribution from the semi-diurnal band, with l(1) = 0.0024, i.e.
   * from latitude dependence. Eq. (9)
   */
  {
    const double shida_l19 = 0.0024e0;
    const double P22 = 3e0 * __cPhi * __cPhi;
    const double __s2lmL = std::sin(2e0 * (lambda - Lambda));
    const double __c2lmL = std::cos(2e0 * (lambda - Lambda));
    dt +=
        -.5e0 * __sphi * __cphi * P22 * (__c2lmL * un + __sphi * __s2lmL * ue);
  }

  /* scale and return */
  const double scale = (GM / GMtb) * std::pow(Re / rtb.norm(), 3) * Re;
  return scale * dr;
}
} /* unnamed namespace */

Eigen::Matrix<double,3,1> displacement() noexcept {
  /* step-1 corrections (time domain) */
  // step1 : Sun 
  // step1 : Moon
  // step1 : out-of-phase (+lat.dependence)
  
  /* step-2 corrections (frequency domain) */
}
