#include "solid_earth_tide.hpp"
#include "geodesy/geodesy.hpp"
#include <cmath>
#include <array>
#include <algorithm>

namespace {
/* @brief Third body (Sun or Moon) Solid earth tide geopotential coefficient
 *        corrections, based on IERS 2010, Anelastic Earth.
 * 
 * Note that the geopotential coefficient corrections computed here (i.e. ΔC 
 * and ΔS) are **added** to the instances provided. That is at output, dC will 
 * be: dC[i](output) = dC[i](input) + dC[i](from this function),
 * and the same goes for the dS array.
 *
 * @param[in] Re Equatorial radius of the Earth [m]
 * @param[in] GM Gravitational constant of Earth
 * @param[in] r_tb ECEF, cartesian position vector of third body (sun or
 *            Moon) [m]
 * @param[in] GM_tb Gravitational constant of Third body
 * @param[out] dC Array where the computed geopotential corrections ΔC are
 *            added. Expected size is 12, in the order:
 *             dC = C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,C43,C44
 *      [indexes] : 0   1   2   3   4   5   6   7   8   9   10  11
 * @param[out] dS Array where the computed geopotential corrections ΔS are
 *            added. Expected size is 12, in the order:
 *            dS = S20,S21,S22,0,S31,S32,S33,0,S41,S42,S43,S44
 * @return Always 0
 */
int iers2010_solid_earth_tide_anelastic_tb(
    double Re, double GM, const Eigen::Matrix<double, 3, 1> &rtb, double GMtb,
    std::array<double, 12> &dC, std::array<double, 12> &dS) noexcept {

  /* get spherical coordinates of third body */
  const auto rtb_spherical = dso::cartesian2spherical(rtb);
  /* trigonometric numbers (of third body) */
  const double __sf = std::sin(rtb_spherical.lat());
  const double __sfsq = __sf * __sf;
  const double __cf = std::cos(rtb_spherical.lat());

  /* compute normalized associated Lagrange polynomials for n=2,3 */
  const double Pnm20 = std::sqrt(5e0) * 0.5e0 * (3e0 * __sfsq - 1e0); // P20
  const double Pnm21 =
      std::sqrt(5.e0 / 3.e0) * 3.e0 * __sf * std::sqrt(1.e0 - __sfsq);   // P21
  const double Pnm22 = std::sqrt(5.e0 / 12.e0) * 3.e0 * (1.e0 - __sfsq); // P22
  const double Pnm30 =
      std::sqrt(7.e0) * 0.5e0 * (5.e0 * __sfsq * __sf - 3.e0 * __sf); // P30
  const double Pnm31 =
      std::sqrt(7.e0 / 6.e0) * 1.5e0 * (5.e0 * __sfsq - 1.e0) * __cf; // P31
  const double Pnm32 =
      std::sqrt(7.e0 / 60.e0) * 15.e0 * __sf * __cf * __cf; // P32
  const double Pnm33 = std::sqrt(7.e0 / 360.e0) * 15.e0 * std::pow(__cf, 3);

  /* IERS 2010, Table 6.3: Nominal values of solid Earth
   * tide external potential Love numbers. Anelastic
   * Earth n  m  Re(k_nm)  Im(k_nm) k_nm^(+)
   * ----------------------------------
   *       2  0  0.30190    0.0
   *       2  1  0.29830   -0.00144
   *       2  2  0.30102   -0.00130
   *       3  0  0.09300
   *       3  1  0.09300
   *       3  2  0.09300
   *       3  3  0.09400
   */
  const double __cl = std::cos(rtb_spherical.lon());
  const double __sl = std::sin(rtb_spherical.lon());
  const double __c2l = std::cos(2e0 * rtb_spherical.lon());
  const double __s2l = std::sin(2e0 * rtb_spherical.lon());
  const double __c3l = std::cos(3e0 * rtb_spherical.lon());
  const double __s3l = std::sin(3e0 * rtb_spherical.lon());

  /* temporary storage; latter on dCt and dSt will be added to dC and dS. */
  std::array<double, 12> dCt = {0e0}, dSt = {0e0};

  /* order n=2, Eq. (6.6) from IERS 2010, ommiting GMj/GM and (Re/r)^3 */
  constexpr const double fac2 = 1e0 / 5e0;
  /* ΔC20 */ dCt[0] = fac2 * Pnm20 * 0.30190e0;
  /* ΔS20 */ dSt[0] = 0e0;
  /* ΔC21 */ dCt[1] = fac2 * Pnm21 * (0.29830e0 * __cl + (-0.00144e0) * __sl);
  /* ΔC21 */ dSt[1] = fac2 * Pnm21 * (0.29830e0 * __sl - (-0.00144e0) * __cl);
  /* ΔC22 */ dCt[2] = fac2 * Pnm22 * (0.30102e0 * __c2l + (-0.00130e0) * __s2l);
  /* ΔC22 */ dSt[2] = fac2 * Pnm22 * (0.30102e0 * __s2l - (-0.00130e0) * __c2l);

  /* order n = 3 Eq. (6.6) from IERS 2010, ommiting GMj/GM and (Re/r)^3.
   * Note that there is no imaginary part of knm for n = 3
   */
  const double fac3 = Re / rtb_spherical.r() / 7e0;
  /* ΔC30 */ dCt[3] = fac3 * Pnm30 * 0.093e0;
  /* ΔS30 */ dSt[3] = 0e0;
  /* ΔC31 */ dCt[4] = fac3 * Pnm31 * 0.093e0 * __cl;
  /* ΔS31 */ dSt[4] = fac3 * Pnm31 * 0.093e0 * __sl;
  /* ΔC32 */ dCt[5] = fac3 * Pnm32 * 0.093e0 * __c2l;
  /* ΔS32 */ dSt[5] = fac3 * Pnm32 * 0.093e0 * __s2l;
  /* ΔC33 */ dCt[6] = fac3 * Pnm33 * 0.094e0 * __c3l;
  /* ΔS33 */ dSt[6] = fac3 * Pnm33 * 0.094e0 * __s3l;

  /* order n = 4 Eq. (6.7) from IERS 2010, ommiting GMj/GM and (Re/r)^3.
   * Coefficients knm(+) are taken from Table 6.3, and are the same for
   * elastic and anelastic case.
   * Note that coefficients knm(+) do not have an imaginary part.
   */
  /* ΔC40 */ dCt[7] = fac2 * Pnm20 * (-0.00089e0) * __cl;
  /* ΔS40 */ dSt[7] = 0e0;
  /* ΔC41 */ dCt[8] = fac2 * Pnm21 * (-0.00080e0) * __cl;
  /* ΔC41 */ dSt[8] = fac2 * Pnm21 * (-0.00080e0) * __sl;
  /* ΔC42 */ dCt[9] = fac2 * Pnm22 * (-0.00057e0) * __c2l;
  /* ΔC42 */ dSt[9] = fac2 * Pnm22 * (-0.00057e0) * __s2l;

  /* scale and add to output arrays */
  const double fac =
      (GMtb / std::pow(rtb_spherical.r(), 3)) * (std::pow(Re, 3) / GM);
  std::transform(dCt.cbegin(), dCt.cend(), dC.cbegin(), dC.begin(),
                 [fac](double dct, double dc) { return dc + dct * fac; });
  std::transform(dSt.cbegin(), dSt.cend(), dS.cbegin(), dS.begin(),
                 [fac](double dst, double ds) { return ds + dst * fac; });

  return 0;
}
} /* unnamed namespace */

int dso::SolidEarthTide::solid_earth_tide_step1(
    const Eigen::Matrix<double, 3, 1> &rMoon,
    const Eigen::Matrix<double, 3, 1> &rSun, std::array<double, 12> &dC,
    std::array<double, 12> &dS) noexcept {

  /* initialize corrections to zero */
  std::fill(dC.begin(), dC.end(), 0e0);
  std::fill(dS.begin(), dS.end(), 0e0);

  /* start with Sun geopotential corrections */
  iers2010_solid_earth_tide_anelastic_tb(cs.Re(), cs.GM(), rSun, mGMsun, dC,
                                         dS);
  /* add Moon */
  iers2010_solid_earth_tide_anelastic_tb(cs.Re(), cs.GM(), rMoon, mGMmoon, dC,
                                         dS);
  /* all done for step 1 */
  return 0;
}
