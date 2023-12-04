#ifndef __DSO_FUNDARG_IERS2010_HPP__
#define __DSO_FUNDARG_IERS2010_HPP__

#include "datetime/calendar.hpp"
#include "geodesy/units.hpp"

namespace dso {

namespace iers2010 {
/* @brief Compute the lunisolar fundamental arguments.
 *
 * This subroutine computes the lunisolar fundamental arguments. The model
 * used is from Simon et al. (1994) as recommended by the IERS Conventions
 * (2010).  Refer to IERS Conventions (2010) Chapter 5 Sections 5.7.1 - 5.7.2
 * (pp. 57 - 59).
 * Though T is strictly TDB, it is usually more convenient to use TT, which
 * makes no significant difference. Julian centuries since J2000 is
 * (JD - 2451545.0)/36525.
 * The expression used is as adopted in IERS Conventions (2010) and is from
 * Simon et al. (1994).
 * All angular arguments are in radians.
 *
 * @param[in]  t     TT, Julian centuries since J2000
 * @param[out] fargs A 5-element array containing the computed fundamental
 *    arguments, in the following order:
 *    fargs[0] -> l  : Mean anomaly of the Moon [rad]
 *    fargs[1] -> lp : Mean anomaly of the Sun [rad]
 *    fargs[2] -> f  : L - OM [rad]
 *    fargs[3] -> d  : Mean elongation of the Moon from the Sun [rad]
 *    fargs[4] -> om : Mean longitude of the ascending node of the Moon [rad]
 * @return Pointer to the first element in fargs
 */
double *fundarg(double tjc, double *fargs) noexcept;

/* @brief Compute the derivatives of lunisolar fundamental arguments 
 *
 * Derivatives are w.r.t. time t (in Julian centuries). All results are in 
 * units of [rad]/[(julian)century]
 *
 * @param[in]  t     TT, Julian centuries since J2000
 * @param[out] fargs A 5-element array containing the computed derivatives of 
 *             fundamental arguments, in the following order:
 *    fargs[0] -> d(l)/dt  in [rad / Julian Centuries] 
 *    fargs[1] -> d(lp)/dt in [rad / Julian Centuries]
 *    fargs[2] -> d(f)/dt  in [rad / Julian Centuries]
 *    fargs[3] -> d(d)/dt  in [rad / Julian Centuries]
 *    fargs[4] -> d(om)/dt in [rad / Julian Centuries]
 * @return Pointer to the first element in fargs
 */
double *fundarg_derivs(double tjc, double *fargs) noexcept;

/* @brief Fundamental argument, IERS Conventions (2003): mean anomaly of
 *       the Moon.
 *
 * Though t is strictly TDB, it is usually more convenient to use TT, which
 * makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean anomaly of the Moon, [rad] in range [0,2π)
 */
#if defined(__clang__)
__attribute__((optnone))
#elif defined(__GNUC__) || defined(__GNUG__)
__attribute__((optimize("O0")))
#endif
inline double fal03(double t) noexcept {
  const double a = dso::sec2rad(std::fmod(
      485868.249036e0 +
          t * (1717915923.2178e0 +
               t * (31.8792e0 + t * (0.051635e0 + t * (-0.00024470e0)))),
      dso::TURNAS));
  return a;
}

/* @brief Fundamental argument, IERS Conventions (2003): mean anomaly of the
 *        Sun.
 *
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean anomaly of the Sun, [rad] in range [0,2π)
 */
inline double falp03(double t) noexcept {
  const double a = dso::sec2rad(std::fmod(
      1287104.793048e0 +
          t * (129596581.0481e0 +
               t * (-0.5532e0 + t * (0.000136e0 + t * (-0.00001149e0)))),
      dso::TURNAS));
  return a;
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of the
 *        Moon minus mean longitude of the ascending node
 *
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return F, [rad] in range [0,2π)
 */
inline double faf03(double t) noexcept {
  const double a = dso::sec2rad(std::fmod(
      335779.526232e0 +
          t * (1739527262.8478e0 +
               t * (-12.7512e0 + t * (-0.001037e0 + t * (0.00000417e0)))),
      dso::TURNAS));
  return a;
}

/* @brief Fundamental argument, IERS Conventions (2003): mean elongation of
 *        the Moon from the Sun.
 *
 * Though t is strictly TDB, it is usually more convenient to use TT, which
 * makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return  D, [rad] in range [0,2π)
 */
inline double fad03(double t) noexcept {
  const double a = dso::sec2rad(std::fmod(
      1072260.703692e0 +
          t * (1602961601.2090e0 +
               t * (-6.3706e0 + t * (0.006593e0 + t * (-0.00003169e0)))),
      dso::TURNAS));
  return a;
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of the
 *        Moon's ascending node.
 *
 * Though t is strictly TDB, it is usually more convenient to use TT, which
 * makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return Ω, [rad] in range [0,2π)
 */
inline double faom03(double t) noexcept {
  const double a = dso::sec2rad(std::fmod(
      450160.398036e0 +
          t * (-6962890.5431e0 +
               t * (7.4722e0 + t * (0.007702e0 + t * (-0.00005939e0)))),
      dso::TURNAS));
  return a;
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Mercury.
 *
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Mercury, [rad] in range [0,2π)
 */
inline double fame03(double t) noexcept {
  return std::fmod(4.402608842e0 + 2608.7903141574e0 * t, dso::D2PI);
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Venus.
 *
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Venus, [rad] in range [0,2π)
 */
inline double fave03(double t) noexcept {
  return std::fmod(3.176146697e0 + 1021.3285546211e0 * t, dso::D2PI);
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Earth.
 *
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Earth, [rad] in range [0,2π)
 */
inline double fae03(double t) noexcept {
  return std::fmod(1.753470314e0 + 628.3075849991e0 * t, dso::D2PI);
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Mars.
 *
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Mars, [rad] in range [0,2π)
 */
inline double fama03(double t) noexcept {
  return std::fmod(6.203480913e0 + 334.0612426700e0 * t, dso::D2PI);
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Jupiter
 *
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Jupiter, [rad] in range [0,2π)
 */
inline double faju03(double t) noexcept {
  return std::fmod(0.599546497e0 + 52.9690962641e0 * t, dso::D2PI);
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Saturn
 *
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Saturn, [rad] in range [0,2π)
 */
inline double fasa03(double t) noexcept {
  return std::fmod(0.874016757e0 + 21.3299104960e0 * t, dso::D2PI);
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Uranus
 *
 * Though t is strictly TDB, it is usually more convenient to use TT,
 * which makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Uranus, [rad] in range [0,2π)
 */
inline double faur03(double t) noexcept {
  return std::fmod(5.481293872e0 + 7.4781598567e0 * t, dso::D2PI);
}

/* @brief Fundamental argument, IERS Conventions (2003): mean longitude of
 *        Neptune
 *
 * Though t is strictly TDB, it is usually more convenient to use TT, which
 * makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return mean longitude of Neptune, [rad] in range [0,2π)
 */
inline double fane03(double t) noexcept {
  return std::fmod(5.311886287e0 + 3.8133035638e0 * t, dso::D2PI);
}

/* @brief Fundamental argument, IERS Conventions (2003): general accumulated
 *        precession in longitude.
 *
 * Though t is strictly TDB, it is usually more convenient to use TT, which
 * makes no significant difference.
 *
 * @param[in] t Julian centuries since J2000.0 [TT] or [TDB]
 * @return general precession in longitude, [rad] in range [0,2π)
 */
inline double fapa03(double t) noexcept {
  return (0.024381750e0 + 0.00000538691e0 * t) * t;
}
//#pragma GCC pop_options

} /* namespace iers2010 */

/* @brief Compute fundamental arguments (overload using dso::MjdEpoch)
 * @see details::fundarg
 */
inline double *fundarg(const dso::MjdEpoch &tt_mjd, double *fargs) noexcept {
  return iers2010::fundarg(tt_mjd.jcenturies_sinceJ2000(), fargs);
}

/* @brief Compute the derivatives of lunisolar fundamental arguments (w.r.t.
 * time t in Julian centuries).
 * All results are in units of [rad]/[(julian)century]
 */
inline double *fundarg_derivs(const dso::MjdEpoch &tt_mjd,
                           double *fargs) noexcept {
  return iers2010::fundarg_derivs(tt_mjd.jcenturies_sinceJ2000(), fargs);
}

/* @brief Overload of fal03 using a MjdEpoch */
inline double fal03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::fal03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of falp03 using a MjdEpoch */
inline double falp03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::falp03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of faf03 using a MjdEpoch */
inline double faf03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::faf03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of fad03 using a MjdEpoch */
inline double fad03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::fad03(tt.jcenturies_sinceJ2000());
}
/* @brief Overload of faom03 using a MjdEpoch */
inline double faom03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::faom03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of fame03 using a MjdEpoch */
inline double fame03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::fame03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of fave03 using a MjdEpoch */
inline double fave03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::fave03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of fae03 using a MjdEpoch */
inline double fae03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::fae03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of fama03 using a MjdEpoch */
inline double fama03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::fama03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of faju03 using a MjdEpoch */
inline double faju03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::faju03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of fasa03 using a MjdEpoch */
inline double fasa03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::fasa03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of faur03 using a MjdEpoch */
inline double faur03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::faur03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of fane03 using a MjdEpoch */
inline double fane03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::fane03(tt.jcenturies_sinceJ2000());
}

/* @brief Overload of fapa03 using a MjdEpoch */
inline double fapa03(const dso::MjdEpoch &tt) noexcept {
  return iers2010::fapa03(tt.jcenturies_sinceJ2000());
}

} /* namespace dso */

#endif
