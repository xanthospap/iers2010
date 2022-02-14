#ifndef __DSO__IERS_1010__
#define __DSO__IERS_1010__

#include "datetime/dtcalendar.hpp"
#include <cmath>
#ifdef DEBUG
#include <cstdio>
#endif

namespace iers2010 {
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
