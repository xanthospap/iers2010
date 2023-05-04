#ifndef __DSO_FUNDARG_IERS2010_HPP__
#define __DSO_FUNDARG_IERS2010_HPP__

#include "datetime/dtcalendar.hpp"
#include "iersc.hpp"

namespace iers2010 {
namespace details {

/* @brief Compute the lunisolar fundamental arguments.
 * @details This subroutine computes the lunisolar fundamental arguments.
 *          The model used is from Simon et al. (1994) as recommended by the
 *          IERS Conventions (2010).  Refer to IERS Conventions (2010)
 *          Chapter 5 Sections 5.7.1 - 5.7.2 (pp. 57 - 59).
 *
 * @param[in]  t     TT, Julian centuries since J2000 (Note 1)
 * @param[out] fargs A 5-element array containing the computed fundamental
 *                   arguments, in the following order:
 *                   fargs[0] -> l  : Mean anomaly of the Moon [rad] (Note 2)
 *                   fargs[1] -> lp : Mean anomaly of the Sun [rad] (Note 2)
 *                   fargs[2] -> f  : L - OM [rad] (Notes 2 and 3)
 *                   fargs[3] -> d  : Mean elongation of the Moon from the Sun
 *                                     [rad] (Note 2)
 *                   fargs[4] -> om : Mean longitude of the ascending node of
 *                                    the Moon [rad] (Note 2)
 *
 * @note
 *       -# Though T is strictly TDB, it is usually more convenient to use
 *          TT, which makes no significant difference. Julian centuries since
 *          J2000 is (JD - 2451545.0)/36525.
 *       -# The expression used is as adopted in IERS Conventions (2010) and
 *          is from Simon et al. (1994).  Arguments are in radians.
 */
void fundarg(double tjc, double *fargs) noexcept;

/* @brief Compute the derivatives of lunisolar fundamental arguments (w.r.t.
 *        time t in Julian centuries). All results are in units of
 *        [rad]/[(julian)century]
 */
void fundarg_derivs(double tjc, double *fargs) noexcept;

} /* namespace details */

/* @brief Compute fundamental arguments (overload using dso::TwoPartDate)
 * @see details::fundarg
 */
inline void fundarg(const dso::TwoPartDate &tt_mjd, double *fargs) noexcept {
  return details::fundarg(tt_mjd.jcenturies_sinceJ2000(), fargs);
}

/* @brief Compute the derivatives of lunisolar fundamental arguments (w.r.t.
 * time t in Julian centuries).
 * All results are in units of [rad]/[(julian)century]
 */
inline void fundarg_derivs(const dso::TwoPartDate &tt_mjd,
                           double *fargs) noexcept {
  return details::fundarg_derivs(tt_mjd.jcenturies_sinceJ2000(), fargs);
}

/* @brief Fundamental argument, IERS Conventions (2003): mean anomaly of
 *       the Moon.
 * Though t is strictly TDB, it is usually more convenient to use
 * TT, which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean anomaly of the Moon, [rad] in range [0,2π)
 */
inline double fal03(double t) noexcept {
  const double a =
      std::fmod(
          485868.249036e0 +
              t * (1717915923.2178e0 +
                   t * (31.8792e0 + t * (0.051635e0 + t * (-0.00024470e0)))),
          iers2010::TURNAS) *
      iers2010::DAS2R;
  return a;
}

/* @brief Overload of fal03 using a TwoPartDate */
inline double fal03(const dso::TwoPartDate &tt) noexcept {
  return fal03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean anomaly of the
 *        Sun.
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean anomaly of the Sun, [rad] in range [0,2π)
 */
inline double falp03(double t) noexcept {
  const double a =
      std::fmod(
          1287104.793048e0 +
              t * (129596581.0481e0 +
                   t * (-0.5532e0 + t * (0.000136e0 + t * (-0.00001149e0)))),
          iers2010::TURNAS) *
      iers2010::DAS2R;
  return a;
}

/* @brief Overload of falp03 using a TwoPartDate */
inline double falp03(const dso::TwoPartDate &tt) noexcept {
  return falp03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of the
 *        Moon minus mean longitude of the ascending node
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return F, [rad] in range [0,2π)
 */
inline double faf03(double t) noexcept {
  const double a =
      std::fmod(
          335779.526232e0 +
              t * (1739527262.8478e0 +
                   t * (-12.7512e0 + t * (-0.001037e0 + t * (0.00000417e0)))),
          iers2010::TURNAS) *
      iers2010::DAS2R;
  return a;
}

/* @brief Overload of faf03 using a TwoPartDate */
inline double faf03(const dso::TwoPartDate &tt) noexcept {
  return faf03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean elongation of
 *        the Moon from the Sun.
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return  D, [rad] in range [0,2π)
 */
inline double fad03(double t) noexcept {
  const double a =
      std::fmod(
          1072260.703692e0 +
              t * (1602961601.2090e0 +
                   t * (-6.3706e0 + t * (0.006593e0 + t * (-0.00003169e0)))),
          iers2010::TURNAS) *
      iers2010::DAS2R;
  return a;
}

/* @brief Overload of fad03 using a TwoPartDate */
inline double fad03(const dso::TwoPartDate &tt) noexcept {
  return fad03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of the
 *        Moon's ascending node.
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return Ω, [rad] in range [0,2π)
 */
inline double faom03(double t) noexcept {
  const double a =
      std::fmod(
          450160.398036e0 +
              t * (-6962890.5431e0 +
                   t * (7.4722e0 + t * (0.007702e0 + t * (-0.00005939e0)))),
          iers2010::TURNAS) *
      iers2010::DAS2R;
  return a;
}

/* @brief Overload of faom03 using a TwoPartDate */
inline double faom03(const dso::TwoPartDate &tt) noexcept {
  return faom03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Mercury.
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Mercury, [rad] in range [0,2π)
 */
inline double fame03(double t) noexcept {
  return std::fmod(4.402608842e0 + 2608.7903141574e0 * t, iers2010::D2PI);
}

/* @brief Overload of fame03 using a TwoPartDate */
inline double fame03(const dso::TwoPartDate &tt) noexcept {
  return fame03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Venus.
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Venus, [rad] in range [0,2π)
 */
inline double fave03(double t) noexcept {
  return std::fmod(3.176146697e0 + 1021.3285546211e0 * t, iers2010::D2PI);
}

/* @brief Overload of fave03 using a TwoPartDate */
inline double fave03(const dso::TwoPartDate &tt) noexcept {
  return fave03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Earth.
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Earth, [rad] in range [0,2π)
 */
inline double fae03(double t) noexcept {
  return std::fmod(1.753470314e0 + 628.3075849991e0 * t, iers2010::D2PI);
}

/* @brief Overload of fae03 using a TwoPartDate */
inline double fae03(const dso::TwoPartDate &tt) noexcept {
  return fae03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Mars.
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Mars, [rad] in range [0,2π)
 */
inline double fama03(double t) noexcept {
  return std::fmod(6.203480913e0 + 334.0612426700e0 * t, iers2010::D2PI);
}

/* @brief Overload of fama03 using a TwoPartDate */
inline double fama03(const dso::TwoPartDate &tt) noexcept {
  return fama03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Jupiter
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Jupiter, [rad] in range [0,2π)
 */
inline double faju03(double t) noexcept {
  return std::fmod(0.599546497e0 + 52.9690962641e0 * t, iers2010::D2PI);
}

/* @brief Overload of faju03 using a TwoPartDate */
inline double faju03(const dso::TwoPartDate &tt) noexcept {
  return faju03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Saturn
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Saturn, [rad] in range [0,2π)
 */
inline double fasa03(double t) noexcept {
  return std::fmod(0.874016757e0 + 21.3299104960e0 * t, iers2010::D2PI);
}

/* @brief Overload of fasa03 using a TwoPartDate */
inline double fasa03(const dso::TwoPartDate &tt) noexcept {
  return fasa03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Uranus
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Uranus, [rad] in range [0,2π)
 */
inline double faur03(double t) noexcept {
  return std::fmod(5.481293872e0 + 7.4781598567e0 * t, iers2010::D2PI);
}

/* @brief Overload of faur03 using a TwoPartDate */
inline double faur03(const dso::TwoPartDate &tt) noexcept {
  return faur03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Neptune
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Neptune, [rad] in range [0,2π)
 */
inline double fane03(double t) noexcept {
  return std::fmod(5.311886287e0 + 3.8133035638e0 * t, iers2010::D2PI);
}

/* @brief Overload of fane03 using a TwoPartDate */
inline double fane03(const dso::TwoPartDate &tt) noexcept {
  return fane03(tt.jcenturies_sinceJ2000());
}

/* @brief Fundamental argument, IERS Conventions (2003): general accumulated
 *        precession in longitude.
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return general precession in longitude, [rad] in range [0,2π)
 */
inline double fapa03(double t) noexcept {
  return (0.024381750e0 + 0.00000538691e0 * t) * t;
}

/* @brief Overload of fapa03 using a TwoPartDate */
inline double fapa03(const dso::TwoPartDate &tt) noexcept {
  return fapa03(tt.jcenturies_sinceJ2000());
}

} /* namespace iers2010 */

#endif
