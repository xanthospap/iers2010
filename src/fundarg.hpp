#ifndef __DSO_FUNDARG_IERS2010_HPP__
#define __DSO_FUNDARG_IERS2010_HPP__

#include "datetime/dtcalendar.hpp"

namespace iers2010 {
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
void fundarg(double tjc, double *fargs) noexcept;

inline void fundarg(dso::TwoPartDate &tt_mjd, double *fargs) noexcept {
  return fundarg(tt_mjd.jcenturies_sinceJ2000(), fargs);
}

}// namespace iers2010

#endif
