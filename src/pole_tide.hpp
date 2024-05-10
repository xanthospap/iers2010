/** @file
 * Functions and structs to handle pole tide, i.e. both pole tide and ocean 
 * pole tide.
 */
#ifndef __DSO_POLE_TIDE_LOADING_HPP__
#define __DSO_POLE_TIDE_LOADING_HPP__

#include "eop.hpp"

namespace dso {

namespace pole_tide_details {  
struct m12Coeffs {double m1, m2;}; 
} /* namespace pole_tide_details */

class poleTide {
  /** Compute m1 and m2 coeffs (IERS 2010, E.q (25).
   *
   * xp, yp are the pole coordinates (IERS 2010, Ch. 5) in [asec].
   * m12Coeffs returned are in [asec].
   */
  static pole_tide_details::m12Coeffs mcoeffs(const MjdEpoch &tt, double xp,
                                              double yp) noexcept {
    /* secular pole in  [milliarcsec] */
    const auto spole = secular_pole(tt);
    printf("%.15f %.15f\n", spole.xs*1e-3, spole.ys*1e-3);
    /* m1, m2 coefficients in [arcsec] */
    return pole_tide_details::m12Coeffs{xp - spole.xs*1e-3,
                                        spole.ys*1e-3 - yp};
  }

public:

  /** @brief Geopotential coefficients corrections (ΔC_12 and ΔS_12) due to 
   *         (solid earth) pole tide.
   *
   * Model implemented here is described in IERS 2010, Ch. 6.4 "Solid Earth 
   * pole tide".
   *
   * @param[in] tt Time of computation (MJD, TT)
   * @param[in] xp Pole X-coordinates at tt in [arcsec]
   * @param[in] yp Pole Y-coordinates at tt in [arcsec]
   * @param[out] dC21 Geopotential coefficient correction for C_21
   * @param[out] dS21 Geopotential coefficient correction for S_21
   */
  static int stokes_coeffs(const MjdEpoch &tt, double xp, double yp,
                           double &dC21, double &dS21) noexcept {
    const auto m12 = mcoeffs(tt, xp, yp);
    const double m1 = m12.m1; /* m1 in [arcsec] */
    const double m2 = m12.m2; /* m2 in [arcsec] */
    dC21 = -1.333e-9 * (m1 + 0.0115e0 * m2);
    dS21 = -1.333e-9 * (m2 - 0.0115e0 * m1);
    printf("m1=%.9f m2=%.9f c21=%.15f s21=%.15f\n", m1, m2, dC21, dS21);
    return 0;
  }
}; /* class poleTide */

}/* namespace dso */

#endif
