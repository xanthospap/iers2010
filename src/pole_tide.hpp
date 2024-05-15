/** @file
 * Functions and structs to handle pole tide, i.e. both pole tide and ocean 
 * pole tide.
 */
#ifndef __DSO_POLE_TIDE_LOADING_HPP__
#define __DSO_POLE_TIDE_LOADING_HPP__

#include "eop.hpp"
#include "coeff_matrix_2d.hpp"

namespace dso {

namespace pole_tide_details {  
struct m12Coeffs {double m1, m2;}; 
  
  /** Compute m1 and m2 coeffs (IERS 2010, E.q (25).
   *
   * TODO it is nor clear what time-scale the input time should be at.
   * xp, yp are the pole coordinates (IERS 2010, Ch. 5) in [arcsec].
   * m12Coeffs returned are in [arcsec].
   *
   * @param[in] t Time of computation.
   * @param[in] xp X-xoordinate of pole in [arcsec].
   * @param[in] yp Y-xoordinate of pole in [arcsec].
   * @return An m12Coeffs instance holding the secolar pole as (xs, ys) in 
   *         [arcsec].
   */
  static m12Coeffs mcoeffs(const MjdEpoch &t, double xp,
                                              double yp) noexcept {
    /* secular pole in  [milliarcsec] */
    const auto spole = secular_pole(t);
    //printf("%.15f %.15f\n", spole.xs*1e-3, spole.ys*1e-3);
    /* m1, m2 coefficients in [arcsec] */
    return pole_tide_details::m12Coeffs{xp - spole.xs*1e-3,
                                        spole.ys*1e-3 - yp};
  }
} /* namespace pole_tide_details */

class poleTide {
public:

  /** @brief Geopotential coefficients corrections (ΔC_12 and ΔS_12) due to 
   *         (solid earth) pole tide.
   *
   * Model implemented here is described in IERS 2010, Ch. 6.4 "Solid Earth 
   * pole tide".
   *
   * @param[in] t Time of computation
   * @param[in] xp Pole X-coordinates at t in [arcsec]
   * @param[in] yp Pole Y-coordinates at t in [arcsec]
   * @param[out] dC21 Geopotential coefficient correction for C_21
   * @param[out] dS21 Geopotential coefficient correction for S_21
   */
  static int stokes_coeffs(const MjdEpoch &t, double xp, double yp,
                           double &dC21, double &dS21) noexcept {
    const auto m12 = pole_tide_details::mcoeffs(t, xp, yp);
    const double m1 = m12.m1; /* m1 in [arcsec] */
    const double m2 = m12.m2; /* m2 in [arcsec] */
    dC21 = -1.333e-9 * (m1 + 0.0115e0 * m2);
    dS21 = -1.333e-9 * (m2 - 0.0115e0 * m1);
    return 0;
  }
}; /* class poleTide */

class oceanPoleTide {
  constexpr static int MAX_DEGREE_DESAI_2002 = 360;
  constexpr static int MAX_ORDER_DESAI_2002 = 360;
public:
  /** @brief Geopotential coefficients corrections (ΔC_12 and ΔS_12) due to 
   *         ocean pole tide.
   *
   * Model implemented here is described in IERS 2010, Ch. 6.5 "Ocean 
   * pole tide". This function will only compute corrections to C21 and S21 
   * coefficients; according to IERS 2010, these account for approximately 90% 
   * of the ocean pole tide potential.
   *
   * @param[in] t Time of computation
   * @param[in] xp Pole X-coordinates at t in [arcsec]
   * @param[in] yp Pole Y-coordinates at t in [arcsec]
   * @param[out] dC21 Geopotential coefficient correction for C_21
   * @param[out] dS21 Geopotential coefficient correction for S_21
   */
  static int stokes_coeffs(const MjdEpoch &t, double xp, double yp,
                           double &dC21, double &dS21) noexcept {
    const auto m12 = pole_tide_details::mcoeffs(t, xp, yp);
    const double m1 = m12.m1; /* m1 in [arcsec] */
    const double m2 = m12.m2; /* m2 in [arcsec] */
    dC21 = -2.1778e-10 * (m1 - 0.01724e0 * m2);
    dS21 = -1.7232e-10*(m2-0.03365*m1);
    return 0;
  }

  //int stokes_coeffs();

  //int parse_desai02_coeffs(
  //    const char *fn, int max_degree, int max_order,
  //    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &A_real,
  //    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &A_imag,
  //    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &B_real,
  //    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &B_imag) noexcept;
  //int parse_desai02_coeffs(
  //    const char *fn,
  //    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &A_real,
  //    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &A_imag,
  //    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &B_real,
  //    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &B_imag) noexcept 
  //    {
  //  return parse_desai02_coeffs(fn, MAX_DEGREE_DESAI_2002, MAX_ORDER_DESAI_2002,
  //                              A_real, A_imag, B_real, B_imag);
  //    }
}; /* class oceanPoleTide */

}/* namespace dso */

#endif
