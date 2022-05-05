#ifndef __DSO__IERS_1010__
#define __DSO__IERS_1010__

#include "datetime/dtcalendar.hpp"
#include "iersc.hpp"
#include "matvec.hpp"
#include <cmath>
#ifdef DEBUG
#include <cstdio>
#endif

namespace iers2010 {
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

/// @brief Compute the lunisolar fundamental arguments.
/// @details This subroutine computes the lunisolar fundamental arguments.
///          The model used is from Simon et al. (1994) as recommended by the
///          IERS Conventions (2010).  Refer to IERS Conventions (2010)
///          Chapter 5 Sections 5.7.1 - 5.7.2 (pp. 57 - 59).
///          This function is a translation/wrapper for the fortran FUNDARG
///          subroutine, found here :
/// https://iers-conventions.obspm.fr/content/chapter5/software/FUNDARG.F
///
/// @param[in]  t     TT, Julian centuries since J2000 (Note 1)
/// @param[out] fargs A 5-element array containing the computed fundamental
///                   arguments, in the following order:
///                   fargs[0] -> l  : Mean anomaly of the Moon [rad] (Note 2)
///                   fargs[1] -> lp : Mean anomaly of the Sun [rad] (Note 2)
///                   fargs[2] -> f  : L - OM [rad] (Notes 2 and 3)
///                   fargs[3] -> d  : Mean elongation of the Moon from the Sun
///                                     [rad] (Note 2)
///                   fargs[4] -> om : Mean longitude of the ascending node of
///                                    the Moon [rad] (Note 2)
/// @return           An integer value always 0.
///
/// @note
///       -# Though T is strictly TDB, it is usually more convenient to use
///          TT, which makes no significant difference. Julian centuries since
///          J2000 is (JD - 2451545.0)/36525.
///       -# The expression used is as adopted in IERS Conventions (2010) and
///          is from Simon et al. (1994).  Arguments are in radians.
///       -# L in this instance is the Mean Longitude of the Moon. OM is the
///          Mean longitude of the ascending node of the Moon.
///
/// @version 25.02.2010
///
/// @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS
///       Technical Note No. 36, BKG (2010); Chapter 5.5.5
int fundarg(double tjc, double *fargs) noexcept;

#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
int fundarg(const dso::datetime<S> &t, double *fargs) noexcept {
  return fundarg(t.jcenturies_sinceJ2000(), fargs);
}

/// @brief Compute the diurnal lunisolar effect on polar motion.
/// @details This function evaluates the model of polar motion for
///          a nonrigid Earth due to tidal gravitation. This polar motion
///          is equivalent to the so-called "subdiurnal nutation." The model
///          is a sum of a first order polynomial and 25 trigonometric terms
///          (15 long periodic and 10 quasi diurnal) with coefficients given
///          in Table 5.1a of the IERS Conventions (2010).
///          This function is a translation/wrapper for the fortran PMSDNUT2
///          subroutine, found here :
///          http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  rmjd Time expressed as Modified Julian date
/// @param[out] dx   The x component of polar motion expressed in
///                  microarcseconds [μas].
/// @param[out] dy   The y component of polar motion expressed in
///                  microarcseconds [μas].
/// @return          An integer value, always 0.
///
/// @warning In the present version this subroutine neglects the linear trend
///          and the long periodic terms of the expansion, for the reasons
///          explained in Section 5.5.1.1 of the IERS Conventions (2010). If
///          the full expansion is needed, set the parameter iband to 0 instead
///          of 1, that is, replace the decleration of iband in the source code.
///
/// @version 13.10.2011
///
/// @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS
///       Technical Note No. 36, BKG (2010); Chapter 5.5.5
int pmsdnut2(double mjd, double &dx, double &dy) noexcept;

/// @brief Compute the diurnal lunisolar effect on polar motion.
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
int pmsdnut2(const dso::datetime<S> &t, double &dx, double &dy) noexcept {
  return pmsdnut2(t.as_mjd(), dx, dy);
}

/// @brief Compute the subdiurnal librations in UT1.
/// @details This function evaluates the model of subdiurnal libration
///          in the axial component of rotation, expressed by UT1 and LOD.
///          This effect is due to the influence of tidal gravitation on
///          the departures of the Earth's mass distribution from the
///          rotational symmetry, expressed by the non-zonal components of
///          geopotential. The amplitudes have been computed for an
///          elastic Earth with liquid core. The adopted truncation level
///          is 0.033 microseconds in UT1 corresponding to the angular
///          displacement of 0.5 microarcseconds or to 0.015 mm at the
///          planet surface. With this truncation level the model contains
///          11 semidiurnal terms. The coefficients of the model are given
///          in Table 5.1b of the IERS Conventions (2010). This function
///          is a translation/wrapper for the fortran UTLIBR subroutine,
///          found here : http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  rmjd Time expressed as Modified Julian date
/// @param[out] dut1 Incremental UT1 in microseconds [μas]
/// @param[out] dlod Incremental LOD in microseconds per day [μas/day]
/// @return          An integer value, always 0.
///
/// @version 23.06.2010
///
/// @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS
///       Technical Note No. 36, BKG (2010); Chapter 5.5.3.3
int utlibr(double mjd, double &dut1, double &dlod) noexcept;

/// @brief Compute the subdiurnal librations in UT1.
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
int utlibr(const dso::datetime<S> &t, double &dut1, double &dlod) noexcept {
  return utlibr(t.as_mjd(), dut1, dlod);
}

/// @details This subroutine computes the effects of the free core nutation.
///          Please note that the table is updated each year (see Note 4).
///          The parameter N needs to be incremented for each additional
///          year in the table.
/// @param[in] mjd   Modified Julian Date, TDB (Note 1)
/// @param[out]  x   CIP offset x component, in microas [μas]
/// @param[out]  y   CIP offset y component, in microas [μas]
/// @param[out] dx   Uncertainty of x component, in microas [μas]
/// @param[out] dy   Uncertainty of y component, in microas [μas]
/// @return          An integer value, always 0.
///
/// @warning Contains data table to be updated each year,
///          see http://syrte.obspm.fr/~lambert/fcn/
/// @note
///      -#  Though the Modified Julian Date (MJD) is strictly TDB, it is
///          usually more convenient to use TT, which makes no significant
///          difference.
///      -#  The updated table is maintained at the website
///          http://syrte.obspm.fr/~lambert/fcn/.
///
/// @version 19.12.2013
/// @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS
///       Technical Note No. 36, BKG (2010); Chapter 5.5.5
int fcnnut(double mjd, double &x, double &y, double &dx, double &dy) noexcept;

/// @brief Compute corrections to the coordinates of the CIP to account for
/// Free Core Nutation.
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
int fcnnut(const dso::datetime<S> &t, double &x, double &y, double &dx,
           double &dy) noexcept {
  return fcnnut(t.as_mjd(), x, y, dx, dy);
}

/// @brief Compute the angular argument which depends on time for 11 tidal
/// argument calculations.
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
dso::Vector3 dehanttideinel_impl(const dso::Vector3 &xsta,
                                 const dso::Vector3 &xsun,
                                 const dso::Vector3 &xmon,
                                 dso::datetime<dso::milliseconds> t) noexcept;

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
/// @param[in]  xsta   Geocentric position of the station (Note 1)
/// @param[in]  xsun   Geocentric position of the Sun (Note 2)
/// @param[in]  xmon   Geocentric position of the Moon (Note 2)
/// @param[in]  t      Datetime UTC (Notes 3 and 4)
/// @return dxtide Displacement vector (Note 5)
///
/// @note
///     -# The station is in ITRF co-rotating frame.  All coordinates,
///        X, Y, and Z, are expressed in meters.
///     -# The position is in Earth Centered Earth Fixed (ECEF) frame.  All
///        coordinates are expressed in meters.
///     -# Time expressed in Coordinated Universal Time (UTC).
///     -#
///     -# The displacement vector is in the geocentric ITRF.  All components
///        are expressed in meters.
///     -# Parameters jc1 and jc2 constitute the date as Julian Centuries in TT
///        time scale. The actual date is given by the addition jc1+jc2.
///        Either jc1 or jc2 can be set to zero.
///     -# Status: Class 1
///     -# This fucnction is part of the package dehanttideinel, see
///        ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/
///
/// @version 19.12.2016
///
/// @cite
///     - Groten, E., 2000, Geodesists Handbook 2000, Part 4,
///     http://www.gfy.ku.dk/~iag/HB2000/part4/groten.htm. See also
///     "Parameters of Common Relevance of Astronomy, Geodesy, and
///     Geodynamics", J. Geod., 74, pp. 134-140
///
///     - Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, "Tidal station
///     displacements", J. Geophys. Res., 102(B9), pp. 20,469-20,477
///
///     - Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
///     IERS Technical Note No. 36, BKG (2010)
///
///     - Pitjeva, E. and Standish, E. M., 2009, "Proposals for the masses
///     of the three largest asteroids, the Moon-Earth mass ratio and the
///     Astronomical Unit", Celest. Mech. Dyn. Astr., 103, pp. 365-372
///
///     - Ries, J. C., Eanes, R. J., Shum, C. K. and Watkins, M. M., 1992,
///     "Progress in the Determination of the Gravitational Coefficient
///     of the Earth", Geophys. Res. Lett., 19(6), pp. 529-531
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
dso::Vector3 dehanttideinel(const dso::Vector3 &xsta, const dso::Vector3 &xsun,
                            const dso::Vector3 &xmon,
                            const dso::datetime<S> &t) noexcept {
  return dehanttideinel_impl(xsta, xsun, xmon,
                             t.template cast_to<dso::milliseconds>());
}

/// Compute tidal corrections of station displacements caused by lunar and
/// solar gravitational attraction.
[[deprecated("use the version with Vector3 instead")]] int
dehanttideinel(const double *, const double *, const double *,
               dso::datetime<dso::seconds>, double *);

/// @brief Compute the diurnal and semi-diurnal variations in Earth Orientation
/// Parameters from ocean tides.
int ortho_eop(double time, double &dx, double &dy, double &dut1);

/// @brief Compute the diurnal and semi-diurnal variations in Earth Orientation
/// Parameters from ocean tides.
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
int ortho_eop(const dso::datetime<S> &t, double &dx, double &dy,
              double &dut1) noexcept {
  return ortho_eop(t.as_mjd(), dx, dy, dut1);
}

/// @brief Evaluate the effects of zonal Earth tides on the rotation of the
///        Earth.
int rg_zont2(double t, double &dut, double &dlod, double &domega) noexcept;
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
int rg_zont2(const dso::datetime<S> &t, double &dut, double &dlod,
             double &domega) noexcept {
  return rg_zont2(t.jcenturies_sinceJ2000(), dut, dlod, domega);
}

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
int cnmtx(double dmjd, double *h) noexcept;

/// @brief Compute the time dependent part of second degree diurnal and
/// semidiurnal tidal potential.
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
int cnmtx(const dso::datetime<S> &t, double *h) noexcept {
  return cnmtx(t.as_mjd(), h);
}

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
/// @param[out] cor_x Tidal correction in x [μas]
/// @param[out] cor_y Tidal correction in y [μas]
int pm_gravi(double t, double &cor_x, double &cor_y) noexcept;

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
/// @param[out] cor_x Tidal correction in x [μas]
/// @param[out] cor_y Tidal correction in y [μas]
/// @param[out] cor_ut1 Tidal correction in UT1-UTC [μsec]
/// @param[out] cor_lod Tidal correction in length of day [μsec]
int pmut1_oceans(double tjc, double &cor_x, double &cor_y, double &cor_ut1,
                 double &cor_lod) noexcept;

/// @brief Perform Lagrangian interpolation.
/// This function performs lagrangian interpolation within a set of (x,y)
/// pairs to give the y value corresponding to xint. This program uses a
/// window of 4 data points to perform the interpolation.
/// For more information see IERS Conventions 2010, sec. 5.5.1.1 "Account of
/// ocean tidal and libration effects in pole coordinates"
/// Translated from the FORTRAN subroutine LAGINT in interp.f found at
/// https://hpiers.obspm.fr/iers/models/
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
/// @return An integer denoting:
///          0 : success
///         -1 : success but xint too close to starting x
///         -2 : success but xint too close to ending x
///         >0 : error
int lagint(const double *x, const double *y, int n, double xint, double &yout,
           int &idx) noexcept;
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
/// @param[in] mjd array of the epochs of data (given in mjd TT) of size n
/// @param[in] x array of x polar motion (millioarcsec) of size n [mas]
/// @param[in] y array of y polar motion (milliarcsec) of size n [mas]
/// @param[in] ut1 array of UT1-UTC (millisec) of size n [ms]
/// @param[in] n number of points in arrays (mjd, x, y, ut1)
/// @param[in] rjd epoch for the interpolated value (given in mjd TT)
/// @param[out] xint interpolated value of x [mas]
/// @param[out] yint interpolated value of y [mas]
/// @param[out] ut1int interpolated value of ut1-utc [ms]
int interp_pole(const double *mjd, const double *x, const double *y,
                const double *ut1, int n, double rjd, double &xint,
                double &yint, double &ut1int) noexcept;

} // namespace iers2010

#endif
