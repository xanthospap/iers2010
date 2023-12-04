#ifndef __DSO_IAU_MODELS_CPP_IMPLEMENTATION_HPP__
#define __DSO_IAU_MODELS_CPP_IMPLEMENTATION_HPP__

#include "datetime/calendar.hpp"
#include "geodesy/units.hpp"

namespace dso {
/* @brief Compute ERA, i.e. Earth rotation angle (IAU 2000 model).
 *
 * The implementation here follows Equation 5.15 in IERS Conventions 2010.
 *
 * @param[in] tt MjdEpoch in [UT1]
 * @return Earth rotation angle [rad], range 0-2pi
 */
double era00(const MjdEpoch &ut1) noexcept;

/* @brief Compute TIO locator s'
 *
 * The TIO locator s', positioning the Terrestrial Intermediate Origin
 * on the equator of the Celestial Intermediate Pole.
 * The TIO locator s' is obtained from polar motion observations by
 * numerical integration, and so is in essence unpredictable. However, it
 * is dominated by a secular drift of about 47 microarcseconds per century,
 * which is the approximation evaluated by the present function.
 *
 * See IERS Conventions 2010, Section 5.5.2
 *
 * @param[in] tt MjdEpoch in [TT]
 * @return the TIO locator s' in [rad]
 */
inline double sp00(const MjdEpoch &tt) noexcept {
  /* Interval between fundamental epoch J2000.0 and current date */
  const double t = tt.jcenturies_sinceJ2000();
  /* Approximate s' */
  return sec2rad(-47e-6 * t);
}

/* @brief X,Y coordinates of the CIP from series based on IAU 2006 precession
 *  and IAU 2000A nutation.
 *
 * The X,Y coordinates are those of the unit vector towards the celestial
 * intermediate pole, in the GCRS. They represent the combined effects of
 * frame bias, precession and nutation.
 * This routine is used for the so called 'CIO-based' transformation, see
 * IERS 2010, 5.5.4. The model used here, uses the full development for the X
 * and Y components and is valid at the μarcsec level (IERS 2010, Section
 * 5.5.4).
 *
 * @param[in] tt MjdEpoch in [TT]
 * @param[out] x CIP X coordinate [rad]
 * @param[out] y CIP Y coordinate [rad]
 * @param[out] fargs If not NULL, at output it will hold the luni-solar and
 *             planetary arguments used in the computation of (X,Y). Since we
 *             are computing them, we might as well return them! If not NULL,
 *             the array should be large enough to hold 14 doubles, i.e.
 *         [l, l', F, D, Om, L_Me, L_Ve, L_E, L_Ma, L_J, L_Sa, L_U, L_Ne, p_A]
 *             all in units of [rad].
 */
int xycip06a(const dso::MjdEpoch &tt, double &xcip, double &ycip,
             double *fargs = nullptr) noexcept;

/** CIO locator compatible with the IAU 2006/2000A precession-nutation model.
 *
 * We follow here IERS 2010 conventions, adopting the development for s(t),
 * or rather the quantity s + XY/2, in the last paragraph of Section 5.5.6.
 * The model is thoroughly described in Table 5.2c --here we do not use the
 * version described in table 5.2d--; the first offers better precision.
 *
 * @param[in] tt MjdEpoch in [TT]
 * @param[in] xcip CIP X coordinate [rad]
 * @param[in] ycip CIP Y coordinate [rad]
 * @return CIO in [rad]
 */
double s06(const dso::MjdEpoch &tt, double xcip, double ycip) noexcept;

/** CIO locator compatible with the IAU 2006/2000A precession-nutation model.
 *
 * We follow here IERS 2010 conventions, adopting the development for s(t),
 * or rather the quantity s + XY/2, in the last paragraph of Section 5.5.6.
 * The model is thoroughly described in Table 5.2c --here we do not use the
 * version described in table 5.2d--; the first offers better precision.
 *
 * @param[in] tt MjdEpoch in [TT]
 * @param[in] xcip CIP X coordinate [rad]
 * @param[in] ycip CIP Y coordinate [rad]
 * @param[in] fargs14 Holds the luni-solar and planetary arguments used in the
 *             computation of CIO (and Xcip, Ycip). It often happens that
 *             these quantities are already computed (e.g. from a call to
 *             xycip06a), so instead of being re-computed within the function,
 *             they can be provided as (input) parameters.
 *             The array is expected to hold the 14 doubles:
 *         [l, l', F, D, Om, L_Me, L_Ve, L_E, L_Ma, L_J, L_Sa, L_U, L_Ne, p_A]
 *             all in units of [rad].
 * @return CIO in [rad]
 */
double s06(const dso::MjdEpoch &tt, double xcip, double ycip,
           const double *const fargs14) noexcept;

namespace detail {
/** @brief X CIP coordinate in GCRS in [μas].
 *
 * The computation follows the IAU 2006 precession and IAU 2000A_R06
 * nutation model(s).
 * Development is valid at the microarcsecond level, based on the IAU 2006
 * precession and IAU 2000A nutation. For more information, see IERS 2010,
 * Section 5.5.4.
 * See also the Table 5.2a published by IERS (2010), i.e.
 * https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2a.txt
 *
 * @param[in] fargs Luni-solar and planetary arguments in the order:
 *  [l, l', F, D, Om, L_Me, L_Ve, L_E, L_Ma, L_J, L_Sa, L_U, L_Ne, p_A]
 *  (i.e. size=14). All units are [rad]
 * @param[in] t = (TT − 2000 January 1d 12h TT) in days/3652 (see IERS 2010,
 *  Sec. 5.3 Eq. 5.2
 * @return CIP X coordinate in GCRS in [μas]
 */
double xcip06a(const double *const fargs, double t) noexcept;

/** @brief Y CIP coordinate in GCRS in [μas].
 *
 * The computation follows the IAU 2006 precession and IAU 2000A_R06
 * nutation model(s).
 * Development is valid at the microarcsecond level, based on the IAU 2006
 * precession and IAU 2000A nutation. For more information, see IERS 2010,
 * Section 5.5.4.
 * See also the Table 5.2b published by IERS (2010), i.e.
 * https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2b.txt
 *
 * @param[in] fargs Luni-solar and planetary arguments in the order:
 *  [l, l', F, D, Om, L_Me, L_Ve, L_E, L_Ma, L_J, L_Sa, L_U, L_Ne, p_A]
 *  (i.e. size=14). All units are [rad]
 * @param[in] t = (TT − 2000 January 1d 12h TT) in days/3652 (see IERS 2010,
 *  Sec. 5.3 Eq. 5.2
 * @return CIP Y coordinate in GCRS in [μas]
 */
double ycip06a(const double *const fargs, double t) noexcept;

/** CIO locator compatible with the IAU 2006/2000A precession-nutation model.
 *
 * We follow here IERS 2010 conventions, adopting the development for s(t),
 * or rather the quantity s + XY/2, in the last paragraph of Section 5.5.6.
 * The model is thoroughly described in Table 5.2c --here we do not use the
 * version described in table 5.2d--; the first offers better precision.
 * In essence, we compute the quantity s + XY/2 using the development 
 * described in Table 5.2c and subsequently subtract the XY/2 (where X is 
 * \p xcip and Y is \p ycip); this is why these quantities are needed.
 *
 * @param[in] fargs Holds the luni-solar and planetary arguments used in the
 *             computation of CIO (and Xcip, Ycip). It often happens that
 *             these quantities are already computed (e.g. from a call to
 *             xycip06a), so instead of being re-computed within the function,
 *             they can be provided as (input) parameters.
 *             The array is expected to hold the 14 doubles:
 *         [l, l', F, D, Om, L_Me, L_Ve, L_E, L_Ma, L_J, L_Sa, L_U, L_Ne, p_A]
 *             all in units of [rad].
 * @param[in] t = (TT − 2000 January 1d 12h TT) in days/3652 (see IERS 2010,
 *  Sec. 5.3 Eq. 5.2
 * @param[in] xcip CIP X coordinate [rad]
 * @param[in] ycip CIP Y coordinate [rad]
 * @return CIO in [rad]
 */
double s06(const double *const fargs, double t, double xcip,
           double ycip) noexcept;
} /* namespace detail */

} /* namespace dso */

#endif
