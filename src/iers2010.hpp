#ifndef __DSO__IERS_1010__
#define __DSO__IERS_1010__

#include "datetime/dtcalendar.hpp"
#include "iersc.hpp"
#include "fundarg.hpp"
#include "geodesy/units.hpp"
#include "eigen3/Eigen/Eigen"
#include <cmath>
#include <vector>
#ifdef DEBUG
#include <cstdio>
#endif

namespace iers2010 {

namespace utils {
/* @brief Compute fundamental arguments and GMST+πι and arrange them in the 
 *        fargs array as: [ GMST+pi, l, lp, f, d, om ] 
 * @warning GMST is computed using the IAU 2006 model (see IERS2010, Sec 5.5.7)
 *        This is different than what most IERS-distributed FORTRAN 
 *        implementations do.
 * @param[in] mjd_tt Date in TT
 * @param[in] dut1   ΔUT1 in [sec]. This value should be looked up using e.g.
 *                   IERS products
 * @param[out] fargs GMST+πι and fundamental arguments [rad], in the order:
 *                   GMST+pi, l, lp, f, d, om
 */
int fargs(const dso::TwoPartDate &mjd_tt, double dut1, double *fargs) noexcept;
} /* iers2010::utils */

dso::TwoPartDate split_fmjd(double fmjd) noexcept;

/// @brief Secular pole coordinates, aka x_s, y_s
/// The coordinates of that secular pole designated (Xs, Ys), given in
/// milliarcseconds. IERS2010, Chapter 7.1.4, Equation (21)
/// @param[in] t Datetime
/// @param[out] xs Secular pole coordinate, X-component in [milliarcseconds]
/// @param[out] ys Secular pole coordinate, Y-component in [milliarcseconds]
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
[[deprecated]]
inline void secular_pole_coordinates(const dso::datetime<S> &t, double &xs,
                                     double &ys) noexcept {
  // t is the date in years of 365.25 days
  constexpr dso::datetime<S> t2000{dso::year(2000), dso::month(1),
                                   dso::day_of_month(1), S(0)};
  constexpr double t2000_mjd = t2000.as_mjd();
  const double dt = (t.as_mjd() - t2000_mjd) / 365.25e0;

  xs = 55e0 + 1.677e0 * dt;
  ys = 320.5e0 + 3.460e0 * dt;
}

/// @brief Effect of Solid earth pole tide.
/// @param[in] t Datetime
/// @param[in] xp Pole coordinates, X-component [arcsec] (see Chapter 5.5.1)
/// @param[in] yp Pole coordinates, Y-component [arcsec] (see Chapter 5.5.1)
/// @param[out] dC21 Changes in the geopotential coefficient C_21 due to the
///                  external potential cused by the solid earth pole tide
/// @param[out] dS21 Changes in the geopotential coefficient S_21 due to the
///                  external potential cused by the solid earth pole tide
/// @see IERS2010, Chapter 6.4
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
[[deprecated]]
void solid_earth_pole_tide(const dso::datetime<S> &t, double xp, double yp,
                           double &dC21, double &dS21) noexcept {
  // get secular pole coordinates -- transform to [arcsec] from milli[arcsec]
  double xs, ys;
  secular_pole_coordinates(t, xs, ys);
  xs *= 1e-3;
  ys *= 1e-3;

  // compute m1 and m2 factors, in seconds of arc
  // Eq. (25), Chapter 7.1.4
  const double m1 = (xp - xs);
  const double m2 = -(yp - ys);

  // equivalent to changes in the geopotential coefficients from tidal
  // deformation (Chapter 6.4)
  dC21 = -1.333e-9 * (m1 + 0.0115 * m2);
  dS21 = -1.333e-9 * (m2 - 0.0115 * m1);
}

/// @brief Effect of Solid earth pole tide.
/// @param[in] t Datetime
/// @param[in] lat Latitude [rad]
/// @param[in] lon Longitude [rad]
/// @param[in] r Geocentric vector (magnitude) directed toward the site [m]
/// @param[in] xp Pole coordinates, X-component [arcsec] (see Chapter 5.5.1)
/// @param[in] yp Pole coordinates, Y-component [arcsec] (see Chapter 5.5.1)
/// @param[out] dC21 Changes in the geopotential coefficient C_21 due to the
///                  external potential cused by the solid earth pole tide
/// @param[out] dS21 Changes in the geopotential coefficient S_21 due to the
///                  external potential cused by the solid earth pole tide
/// @return Perturbation caused by the pole tide
/// @see IERS2010, Chapter 6.4
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
[[deprecated]]
double solid_earth_pole_tide(const dso::datetime<S> &t, double lat, double lon,
                             double r, double xp, double yp, double &dC21,
                             double &dS21) noexcept {
  // get secular pole coordinates -- transform to [arcsec] from milli[arcsec]
  double xs, ys;
  secular_pole_coordinates(t, xs, ys);
  xs *= 1e-3;
  ys *= 1e-3;

  // compute m1 and m2 factors, in seconds of arc
  // Eq. (25), Chapter 7.1.4
  const double m1 = (xp - xs);
  const double m2 = -(yp - ys);

  // equivalent to changes in the geopotential coefficients from tidal
  // deformation (Chapter 6.4)
  dC21 = -1.333e-9 * (m1 + 0.0115 * m2);
  dS21 = -1.333e-9 * (m2 - 0.0115 * m1);

  const double m1c = m1 * std::cos(lon);
  const double m2s = m2 * std::sin(lon);
  const double s2t = std::sin(2e0 * (DPI / 2e0 - lat));
  const double factor = s2t * (m1c + m2s);

  return -(OmegaEarth * OmegaEarth) * (r * r * factor) / 2e0;
}

/* @brief Compute the diurnal lunisolar effect on polar motion.
 *
 * This function evaluates the model of polar motion for a nonrigid Earth due 
 * to tidal gravitation. This polar motion is equivalent to the so-called 
 * "subdiurnal nutation." The model is a sum of a first order polynomial and 
 * 25 trigonometric terms (15 long periodic and 10 quasi diurnal) with 
 * coefficients given in Table 5.1a of the IERS Conventions (2010). This 
 * function is a translation/wrapper for the fortran PMSDNUT2 subroutine, 
 * found here :
 *          http://maia.usno.navy.mil/conv2010/software.html
 *
 * @param[in] mjd_tt Date in TT
 * @param[in] dut1  ΔUT1 in [sec]. This value should be looked up using e.g.
 *                  IERS products
 * @param[out] dx   The x component of polar motion expressed in
 *                  microarcseconds [μas].
 * @param[out] dy   The y component of polar motion expressed in
 *                  microarcseconds [μas].
 * @return          An integer value, always 0.
 *
 * @warning In the present version this subroutine neglects the linear trend
 * and the long periodic terms of the expansion, for the reasons explained in 
 * Section 5.5.1.1 of the IERS Conventions (2010). If the full expansion is 
 * needed, set the parameter iband to 0 instead of 1, that is, replace the 
 * decleration of iband in the source code.
 */
int pmsdnut2(const dso::TwoPartDate &mjd_tt, double dut1, double &dx,
             double &dy) noexcept;

/* @brief Overload of pmsdnut2, given the Fundamental (Delaunay) arguments and
 *        GMST+π at the requested epoch (in fargs).
 *
 * @param[in]  mjd_tt Date in TT
 * @param[out] fargs GMST+πι and fundamental arguments, in the order:
 *             GMST+pi, l, lp, f, d, om (see the function utils::fargs).
 * @param[out] dx   The x component of polar motion expressed in
 *                  microarcseconds [μas].
 * @param[out] dy   The y component of polar motion expressed in
 *                  microarcseconds [μas].
 * @return          An integer value, always 0.
 */
int pmsdnut2(const dso::TwoPartDate &mjd_tt, const double *const fargs,
             double &dx, double &dy) noexcept;

/* @brief Compute the subdiurnal librations in UT1.
 *
 * @details This function evaluates the model of subdiurnal libration
 *    in the axial component of rotation, expressed by UT1 and LOD. This 
 *    effect is due to the influence of tidal gravitation on the departures of 
 *    the Earth's mass distribution from the rotational symmetry, expressed by 
 *    the non-zonal components of geopotential. The amplitudes have been 
 *    computed for an elastic Earth with liquid core. The adopted truncation 
 *    level is 0.033 microseconds in UT1 corresponding to the angular 
 *    displacement of 0.5 microarcseconds or to 0.015 mm at the planet surface. 
 *    With this truncation level the model contains 11 semidiurnal terms. The 
 *    coefficients of the model are given in Table 5.1b of the IERS 
 *    Conventions (2010). This function is a translation/wrapper for the 
 *    fortran UTLIBR subroutine, found here: 
 *    http://maia.usno.navy.mil/conv2010/software.html
 *
 *    Note that contrary to the FORTRAN implementation, we expect here an 
 *    extra input argument, i.e. dut1in. That is because in this version, 
 *    GMST is computed using the IAU2006(A) model, which depends on both TT 
 *    and UT1.
 *
 * @param[in] mjd_tt  Date in TT
 * @param[in] dut1in  ΔUT1 in [sec]. This value should be looked up using e.g.
 *                    IERS products
 * @param[out] dut1 Incremental UT1 in microseconds [μas]
 * @param[out] dlod Incremental LOD in microseconds per day [μas/day]
 * @return          An integer value, always 0.
 *
 * @version 23.06.2010
 *
 * @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS
 *       Technical Note No. 36, BKG (2010); Chapter 5.5.3.3
 */
int utlibr(const dso::TwoPartDate &mjd_tt, double dut1in, double &dut1,
           double &dlod) noexcept;

/* @brief Compute the subdiurnal librations in UT1.
 * @details This function evaluates the model of subdiurnal libration
 *          in the axial component of rotation, expressed by UT1 and LOD.
 *          This effect is due to the influence of tidal gravitation on
 *          the departures of the Earth's mass distribution from the
 *          rotational symmetry, expressed by the non-zonal components of
 *          geopotential. The amplitudes have been computed for an
 *          elastic Earth with liquid core. The adopted truncation level
 *          is 0.033 microseconds in UT1 corresponding to the angular
 *          displacement of 0.5 microarcseconds or to 0.015 mm at the
 *          planet surface. With this truncation level the model contains
 *          11 semidiurnal terms. The coefficients of the model are given
 *          in Table 5.1b of the IERS Conventions (2010). This function
 *          is a translation/wrapper for the fortran UTLIBR subroutine,
 *          found here : http://maia.usno.navy.mil/conv2010/software.html
 *
 * @param[in]  mjd_tt Date in TT
 * @param[out] fargs GMST+πι and fundamental arguments, in the order:
 *             GMST+pi, l, lp, f, d, om (see the function utils::fargs).
 * @param[out] dut1 Incremental UT1 in microseconds [μas]
 * @param[out] dlod Incremental LOD in microseconds per day [μas/day]
 * @return          An integer value, always 0.
 *
 * @version 23.06.2010
 *
 * @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS
 *       Technical Note No. 36, BKG (2010); Chapter 5.5.3.3
 */
int utlibr(const dso::TwoPartDate &mjd_tt, const double *const fargs,
           double &dut1, double &dlod) noexcept;

/* @details This subroutine computes the effects of the free core nutation.
 * @param[in] mjd   Date in TDB; note howver that usually more it is more 
 *                  convenient to use TT, which makes no significant 
 *                  difference.
 * @param[out]  x   CIP offset x component, in microas [μas]
 * @param[out]  y   CIP offset y component, in microas [μas]
 * @param[out] dx   Uncertainty of x component, in microas [μas]
 * @param[out] dy   Uncertainty of y component, in microas [μas]
 * @return          An integer value, always 0.
 *
 * @warning Contains data table to be updated each year,
 *          see http://syrte.obspm.fr/~lambert/fcn/
 *
 * @version 19.12.2013
 * @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS
 *       Technical Note No. 36, BKG (2010); Chapter 5.5.5
 */
int fcnnut(const dso::TwoPartDate &mjd_tt, double &x, double &y, double &dx,
           double &dy) noexcept;

/// @brief Compute the angular argument which depends on time for 11 tidal
/// argument calculations.
/// @details The purpose of the function is to compute the angular astronomical
///          argument, which depends on time, for 11 tidal argument
///          calculations. The order of the 11 angular quantities in vector
///          angle are given below: 01-M2, 02-S2, 03-N2, 04-K2, 05-K1, 06-O1,
///          07-P1, 08-Q1, 09-Mf, 10-Mm, 11-Ssa (See Reference 1) This function
///          is a translation/wrapper for the fortran ARG2 subroutine, found
///          here : http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  iyear Four digit year (Note 1)
/// @param[in]  day   (Fractional) Day of Year Greenwich Time (Note 2)
/// @param[out] angle Angular argument for Schwiderski computation, in [rad]
///                   (Notes 3, 4 and 5). Array of size > 11.
/// @return           An integer value which can be:
///                   Returned Value | Status
///                   ---------------|--------------------
///                               -1 | Error; Invalid year
///                                0 | All ok
///
/// @note
///  -# This subroutine is valid only after 1973 CE. 
///  -# Example: 32.5 for February 1 12 Noon
///     Example: 1.25 for January 1 6 AM
///  -# Ocean loading phases computed from Schwiderski's models
///     refer to the phase of the associated solid Earth tide
///     generating potential at the zero meridian according to <br>
///      OL_DR = OL_AMP ' COS (SE_PHASE" - OL_PHASE) <br>
///     where OL = OCEAN LOADING TIDE,<br>
///           SE = SOLID EARTH TIDE GENERATING POTENTIAL.<br>
///     If the harmonic tide development of Cartwright, et al.
///     (CTE) (1971, 1973) is used, make sure that SE_PHASE"
///     take into account:
///     - the sign of SE_AMP in the tables of Cartwright et al.
///     - that CTE'S SE_PHASE refers to a sine rather than a
///       cosine function if (N+M) = (DEGREE + ORDER) of the tide
///       spherical harmonic is odd.
///     i.e. SE_PHASE" = TAU(T) ' N1 + S(T) ' N2 + H(T) ' N3 <br>
///                 + P(T) ' N4 + N'(T) ' N5 + PS(T) ' N6 <br>
///                 + PI   If CTE'S amplitude coefficient < 0 <br>
///                 + PI/2 If (DEGREE + N1) is odd <br>
///     where TAU ... PS = astronomical arguments,<br>
///           N1 ... N6 = CTE'S argument numbers.<br>
///     Most tide generating software compute SE_PHASE" (for use
///     with cosines).
///  -# The input array angle must be able to hold at least 11 doubles.
///
/// @version  07.10.2011
///
/// @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS
///       Technical Note No. 36, BKG (2010)
/// @cite Schwiderski, E., 1983, "Atlas of Ocean Tidal Charts and Maps, Part I:
///       The Semidiurnal Principal Lunar Tide M2," Marine Geodesy, 6,
///       pp. 219-256.
int arg2(int iyear, double fdoy, double *angles) noexcept;

/// @brief Compute the angular argument which depends on time for 11 tidal
/// argument calculations.
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
int arg2(const dso::datetime<S> &t, double *angles) noexcept {
  const dso::ydoy_date ydoy = t.as_ydoy();
  const double day_fraction = t.sec().fractional_days();
  return arg2(ydoy.__year.as_underlying_type(),
              day_fraction +
                  static_cast<double>(ydoy.__doy.as_underlying_type()),
              angles);
}

/// @brief Compute tidal corrections of station displacements caused by lunar
/// and solar gravitational attraction. This is just the implementation, use the
/// generic template function instead.
int dehanttideinel_impl(
    double julian_centuries_tt, double fhr_ut,
    const Eigen::Matrix<double, 3, 1> &xsun,
    const Eigen::Matrix<double, 3, 1> &xmon,
    const std::vector<Eigen::Matrix<double, 3, 1>> &xsta_vec,
    std::vector<Eigen::Matrix<double, 3, 1>> &xcor_vec) noexcept;

/// @details This function computes the station tidal displacement
///          caused by lunar and solar gravitational attraction (see
///          References). The computations are calculated by the following
///          steps:<br> <b>Step 1):</b> General degree 2 and degree 3
///          corrections
///          + CALL ST1IDIU + CALL ST1ISEM + CALL ST1L1.<br>
///          <b>Step 2):</b> CALL STEP2DIU + CALL STEP2LON<br>
///          It has been decided that the <b>Step 3</b> non-correction for
///          permanent tide would not be applied in order to avoid a jump in the
///          reference frame. This Step 3 must be added in order to get the
///          non-tidal station position and to conform with the IAG Resolution.
///          Details for the implementation can be found in IERS2010, Chapter
///          "7.1.1 Effects of the solid Earth tides".
///
///          This function is a translation/wrapper for the fortran
///          DEHANTTIDEINEL subroutine, found here :
///          http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  t     Datetime UTC
/// @param[in]  xsta_vec  A vector of ECEF positions of the input stations,
///                    (X,Y,Z) in [m]
/// @param[in]  xsun  ECEF position of the Sun (X,Y,Z) in [m]
/// @param[in]  xmon  ECEF position of the Moon (X,Y,Z) in [m]
/// @param[out] xcor_vec Displacements for the sites passed in through the 
///                   xsta_vec vector (in the same order). Corrections are in
///                   the (input) ECEF, as (dX, dY, dZ) in [m]
/// @return     Always 0
///
/// @note This fucnction is part of the package dehanttideinel, see
///   ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/
///
/// @version 19.12.2016
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
int
dehanttideinel(
    const dso::datetime<S>& tutc,
    const Eigen::Matrix<double, 3, 1> &xsun,
    const Eigen::Matrix<double, 3, 1> &xmon,
    const std::vector<Eigen::Matrix<double, 3, 1>> &xsta_vec,
    std::vector<Eigen::Matrix<double, 3, 1>> &xcor_vec) noexcept
{
  // copy date, going to translate UTC to TT
  auto t = tutc;
  
  // need to convert UTC to TT: TT = UTC + 32.184[sec] + DAT =
  int dat = dso::dat(t.mjd());
  dso::milliseconds msec(dso::milliseconds(dat * 1e3) +
                         dso::milliseconds(32184));
  t.add_seconds(dso::cast_to<dso::milliseconds, dso::nanoseconds>(msec));

  // fractional hours in UTC day
  const double fhr = tutc.sec().to_fractional_seconds() / 3600e0;

  return dehanttideinel_impl(t.jcenturies_sinceJ2000(), fhr, xsun, xmon,
                             xsta_vec, xcor_vec);
}

/* Never use this, it is just for debugging/testing purposes */
/// Compute tidal corrections of station displacements caused by lunar and
/// solar gravitational attraction.
[[deprecated("use the version with Vector3 instead")]] int
dehanttideinel(const double *, const double *, const double *,
               dso::datetime<dso::seconds>, double *);

/// @brief Compute the diurnal and semi-diurnal variations in EOPs from ocean 
///        tides.
/// @details  The purpose of the function is to compute the diurnal and semi-
///           diurnal variations in Earth Orientation Parameters (x,y, UT1)
///           from ocean tides.
///           This function is a translation/wrapper for the fortran ORTHO_EOP
///           subroutine, found here :
///           http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  dmjd  Modified Julian Date
/// @param[out] dx    delta_x, in microarcseconds [μas]
/// @param[out] dy    delta_y, in microarcseconds [μas]
/// @param[out] dut1  delta_UT1, in microseconds [μs]
/// @return           An integer, always 0.
///
/// @note
///    -# The diurnal and semidiurnal orthoweights fit to the 8 constituents
///       are listed in Reference 1.
///
/// @version 19.03.2010
///
/// @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS
///       Technical Note No. 36, BKG (2010); Chapter 5.5.5
/// @cite Ray, R. D., Steinberg, D. J., Chao, B. F., and Cartwright, D. E.,
///       "Diurnal and Semidiurnal Variations in the Earth's Rotation
///        Rate Induced by Ocean Tides", 1994, Science, 264, pp. 830-832
int ortho_eop(const dso::TwoPartDate &mjd, double &dx, double &dy,
              double &dut1) noexcept;
int ortho_eop(double fmjd, double &dx, double &dy, double &dut1) noexcept;

/* @details This function evaluates the effects of zonal Earth tides on the
 *          rotation of the Earth.  The model used is a combination of Yoder
 *          et al. (1981) elastic body tide, Wahr and Bergen (1986) inelastic
 *          body tide, and Kantha et al. (1998) ocean tide models
 *          as recommended by the IERS Conventions (2010).  Refer to
 *          Chapter 8 pp. 105 - 106.  The latest version of the model is
 *          located at http://tai.bipm.org/iers/convupdt/convupdt_c8.html. This
 *          function is a translation/wrapper for the fortran RG_ZONT2
 *          subroutine, found here :
 *          http://maia.usno.navy.mil/conv2010/software.html
 *
 * @param[out] fundarg are the fundamental arguments, in the order:
 *             [l, lp, f, d, om]
 * @param[out] dut    Effect on UT1 [sec]
 * @param[out] dlod   Effect on excess length of day (LOD) [sec/day]. The
 *                    phrase 'per day' is generally understood, so it has
 *                    been omitted commonly in speech and literature. See:
 *                    Stephenson, F. R., Morrison, L. V., Whitrow, G. J., 1984,
 *                    "Long-Term Changes in the Rotation of the Earth: 700
 *                    B. C. to A. D. 1980 [and Discussion]", Phil. Trans. Roy.
 *                    Soc. of London. Series A, 313, pp. 47 - 70.
 * @param[out] domega Effect on rotational speed [rad/sec]
 * @return            An integer value always 0.
 *
 * @version 20.12.2011
 *
 * @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
 *     IERS Technical Note No. 36, BKG (2010)
 */
int rg_zont2(const double *const fundarg, double &dut, double &dlod,
             double &domega) noexcept;

/* @details This function evaluates the effects of zonal Earth tides on the
 *    rotation of the Earth.  The model used is a combination of Yoder et al. 
 *    (1981) elastic body tide, Wahr and Bergen (1986) inelastic body tide, 
 *    and Kantha et al. (1998) ocean tide models as recommended by the IERS 
 *    Conventions (2010). Refer to Chapter 8 pp. 105 - 106. This function is a t
 *    ranslation/wrapper for the fortran RG_ZONT2 subroutine, found here :
 *    http://maia.usno.navy.mil/conv2010/software.html
 *
 *    Note that contrary to the FORTRAN implementation, we expect here an 
 *    extra input argument, i.e. dut1in. That is because in this version, 
 *    GMST is computed using the IAU2006(A) model, which depends on both TT 
 *    and UT1.
 *
 * @param[in]  mjd_tt Date in TT
 * @param[out] dut    Effect on UT1 [sec]
 * @param[out] dlod   Effect on excess length of day (LOD) [sec/day]. The
 *                    phrase 'per day' is generally understood, so it has
 *                    been omitted commonly in speech and literature. See:
 *                    Stephenson, F. R., Morrison, L. V., Whitrow, G. J., 1984,
 *                    "Long-Term Changes in the Rotation of the Earth: 700
 *                    B. C. to A. D. 1980 [and Discussion]", Phil. Trans. Roy.
 *                    Soc. of London. Series A, 313, pp. 47 - 70.
 * @param[out] domega Effect on rotational speed [rad/sec]
 * @return            An integer value always 0.
 *
 * @version 20.12.2011
 *
 * @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
 *     IERS Technical Note No. 36, BKG (2010)
 */
int rg_zont2(const dso::TwoPartDate &tt_mjd, double &dut, double &dlod,
             double &domega) noexcept;

/// Compute the global total FCULa mapping function.
double fcul_a(double, double, double, double) noexcept;

/// Computes the global total FCULb mapping function.
double fcul_b(double, double, double, double) noexcept;

// Determine the total zenith delay following Mendes and Pavlis, 2004.
int fcul_zd_hpa(double, double, double, double, double, double &, double &,
                double &) noexcept;

/// Determines the total zenith delay following (Mendes and Pavlis, 2004).
int fculzd_hpa(double, double, double, double, double, double &, double &,
               double &) noexcept;

/// Compute the Global Mapping Functions (GMF).
int gmf(double, double, double, double, double, double &, double &) noexcept;

/// Compute the Vienna Mapping Functions 1 (VMF1), to be used with "a"
/// coefficients computed for a given site.
int vmf1(double, double, double, double, double, double &, double &) noexcept;

/// Compute the Vienna Mapping Functions 1 (VMF1), with height corrections,
/// to be used with "a" coefficients computed for a grid.
int vmf1_ht(double, double, double, double, double, double, double &,
            double &) noexcept;

/// Compute the Global Pressure and Temperature (GPT), based on spherical
/// harmonics up to degree and order 9.
int gpt(double, double, double, double, double &, double &, double &) noexcept;

/// Compute the Global Pressure and Temperature 2 model (GPT2),
/// and the "a" coefficients for vmf1_ht.
int gpt2(double, const double *, const double *, const double *, int, int,
         double *, double *, double *, double *, double *, double *, double *,
         const char *ifile = nullptr);
int gpt2(double dmjd, double dlat, double dlon, double hell, int it, double &p,
         double &t, double &dt, double &e, double &ah, double &aw, double &undu,
         const char *ifile = nullptr);

namespace oeop {
/// @brief Compute the time dependent part of second degree diurnal and
/// semidiurnal tidal potential.
int cnmtx(const dso::TwoPartDate &mjd, double *h) noexcept;
int cnmtx(double fmjd, double *h) noexcept;

/// @brief Compute the time dependent part of second degree diurnal and
/// semidiurnal tidal potential.
//#if __cplusplus >= 202002L
//template <gconcepts::is_sec_dt S>
//#else
//template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
//#endif
//int cnmtx(const dso::datetime<S> &t, double *h) noexcept {
//  return cnmtx(t.as_mjd(), h);
//}

} // namespace oeop

namespace interp {

/// @brief Compute diurnal lunisolar effet on polar motion (")
/// This function provides, in time domain, the diurnal lunisolar effet on
/// polar motion (").
/// The fundamental lunisolar arguments are those of Simon et al.
/// These corrections should be added to "average" EOP values to get estimates
/// of the instantaneous values.
/// For more information see IERS Conventions 2010, sec. 5.5.1.1 "Account of
/// ocean tidal and libration effects in pole coordinates"
/// Translated from the FORTRAN subroutine PMUT1_GRAVI in interp.f found at
/// https://hpiers.obspm.fr/iers/models/
/// @param[in] tjc Epoch of interest given in Julian Centuries TT
/// @param[out] cor_x Tidal correction in x ['', arcssec]
/// @param[out] cor_y Tidal correction in y ['', arcsex]
int pm_gravi(double tjc, double &cor_x, double &cor_y) noexcept;

/// @brief Compute tidal effets on polar motion.
/// This function provides, in time domain, the diurnal/subdiurnal
/// tidal effets on polar motion ("), UT1 (s) and LOD (s). The tidal terms,
/// listed in the program above, have been extracted from the procedure
/// ortho_eop.f coed by Eanes in 1997.
/// For more information see IERS Conventions 2010, sec. 5.5.1.1 "Account of
/// ocean tidal and libration effects in pole coordinates"
/// Translated from the FORTRAN subroutine PMUT1_OCEANS in interp.f found at
/// https://hpiers.obspm.fr/iers/models/
/// @param[in] tjc Epoch of interest given in Julian Centuries TT
/// @param[out] cor_x Tidal correction in x ['', arcsec]
/// @param[out] cor_y Tidal correction in y ['', arcsec]
/// @param[out] cor_ut1 Tidal correction in UT1-UTC [sec]
/// @param[out] cor_lod Tidal correction in length of day [sec]
int pmut1_oceans(double tjc, double &cor_x, double &cor_y, double &cor_ut1,
                 double &cor_lod) noexcept;

/// @brief Perform Lagrangian interpolation.
///
/// This function performs lagrangian interpolation within a set of (x,y)
/// pairs to give the y value corresponding to xint. This program uses a
/// window of 4 data points to perform the interpolation, but this can change
/// via the order parameter.
///
/// For more information see IERS Conventions 2010, sec. 5.5.1.1 "Account of
/// ocean tidal and libration effects in pole coordinates"
/// Translated from the FORTRAN subroutine LAGINT in interp.f found at
/// https://hpiers.obspm.fr/iers/models/
/// 
/// @param[in] x array of values of the independent variable; size n
/// @param[in] y array of function values corresponding to x; size n
/// @param[in] n number of points, aka size of x and y arrays
/// @param[in] xint the x-value for which estimate of y is desired
/// @param[out] yout the y value returned to caller
/// @param[inout] idx If idx>=0, then this value will be used as the index for
///               the interpolation, aka it will be considered that:
///               xint >= x[idx] && xint < x[idx + 1]                      (1)
///               holds (so that the function can skip the search).
///               If idx<0, then it will be ignored and the function will
///               search for the index.
///               At output, idx will hold the value/index for which (1)
///               holds (so that it can be used later if needed).
/// @param[in] order Lagrangian interpolation order. The number of data points 
///               to be used for the interpolation will be order+1. This 
///               parameter should have an odd, integer value.
/// @return An integer denoting:
///          0 : success
///         >0 : error
int lagint(const double *x, const double *y, int n, double xint, double &yout,
           int &idx, int order = 3) noexcept;
} // namespace interp

/// @brief Account of ocean tidal and libration effects in pole coordinates
/// The subdaily variations are not part of the polar motion values reported
/// to and distributed by the IERS and are therefore to be added after
/// interpolation. This is appropriately done by this function, which
/// interpolates series of X_iers, Y_iers values to a chosen date and then adds
/// the contribution for this date of (i) the tidal terms and (ii) the diurnal
/// components.
/// This function takes a series of x, y, and UT1-UTC values and interpolates
/// them to an epoch of choice. This routine assumes that the values of x and
/// y are in seconds of arc and that UT1-UTC is in seconds of time. At least
/// one point before and one point after the epoch of the interpolation point
/// are necessary in order for the interpolation scheme to work.
/// For more information see IERS Conventions 2010, sec. 5.5.1.1 "Account of
/// ocean tidal and libration effects in pole coordinates"
/// Translated from the FORTRAN subroutine INTERP in interp.f found at
/// https://hpiers.obspm.fr/iers/models/
/// @note The original FORTRAN routine, uses arcseconds/seconds as input/output
///       units. Here we use [mas] and [ms] respectively (the same units used
///       in the Bulletin B/C publications).
/// @param[in] mjd array of the epochs of data as MJD of size n (see note)
/// @param[in] x array of x polar motion (arcsec) of size n [arcsec]
/// @param[in] y array of y polar motion (arcsec) of size n [arcsec]
/// @param[in] ut1 array of UT1-UTC (sec) of size n [sec]
/// @param[in] n number of points in arrays (mjd, x, y, dut)
/// @param[in] rjd epoch for the interpolated value, as MJD (see note)
/// @param[out] xint interpolated value of x [arcsec]
/// @param[out] yint interpolated value of y [arcsec]
/// @param[out] ut1int interpolated value of ut1-utc [sec]
/// @param[out] corlod tidal correction in length of day (LOD) in [sec]
///
/// @note The time scales in params mjd and rjd should be the same.
int interp_pole(const double *mjd, const double *x, const double *y,
                const double *ut1, int n, double rjd, double &xint,
                double &yint, double &ut1int, double &corlod) noexcept;

Eigen::Matrix<double, 3, 1>
relativistic_correction(const Eigen::Matrix<double, 3, 1> &rgcrs,
                        const Eigen::Matrix<double, 3, 1> &vgcrs,
                        double c = iers2010::C, double J = 9.8e8,
                        double gme = iers2010::GMe) noexcept;
Eigen::Matrix<double, 3, 1>
relativistic_correction(const Eigen::Matrix<double, 3, 1> &r,
                        const Eigen::Matrix<double, 3, 1> &Res,
                        const Eigen::Matrix<double, 3, 1> &dResdt,
                        const Eigen::Matrix<double, 3, 1> &v,
                        double c = iers2010::C, double J = 9.8e8,
                        double gme = iers2010::GMe,
                        double gms = iers2010::GMSun) noexcept;
} // namespace iers2010

#endif
