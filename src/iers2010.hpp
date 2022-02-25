#ifndef __DSO__IERS_1010__
#define __DSO__IERS_1010__

#include "datetime/dtcalendar.hpp"
#include "iersc.hpp"
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
int fundarg(double tjc, double *fargs) noexcept;

/// @brief Compute the diurnal lunisolar effect on polar motion.
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

/// @brief Compute corrections to the coordinates of the CIP to account for
/// Free Core Nutation.
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

/// Compute tidal corrections of station displacements caused by lunar and
/// solar gravitational attraction.
int dehanttideinel(const double *, const double *, const double *,
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

} // namespace iers2010

#endif
