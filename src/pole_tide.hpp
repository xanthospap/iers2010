/** @file
 * Functions and structs to handle pole tide, i.e. both pole tide and ocean 
 * pole tide.
 */
#ifndef __DSO_POLE_TIDE_LOADING_HPP__
#define __DSO_POLE_TIDE_LOADING_HPP__

#include "iersconst.hpp"
#include "eop.hpp"
#include "stokes_coefficients.hpp"
#include "geodesy/transformations.hpp"

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
inline m12Coeffs mcoeffs(const MjdEpoch &t, double xp, double yp) noexcept {
  /* secular pole in  [mas] */
  const auto spole = secular_pole(t);
  /* m1, m2 coefficients in [arcsec] */
  return m12Coeffs{xp - spole.xs * 1e-3, spole.ys * 1e-3 - yp};
}

/* max degree of coefficients for the Desai 2002 model. Applies to
 * A_real, A_imag, B_real and B_imag when the degree and order of the
 * instance are > 2.
 */
constexpr int MAX_DEGREE_DESAI_2002 = 360;

/* max order of coefficients for the Desai 2002 model. Applies to
 * A_real, A_imag, B_real and B_imag when the degree and order of the
 * instance are > 2.
 */
constexpr int MAX_ORDER_DESAI_2002 = 360;

/** Parse coefficients of the Desai 2002 model for Ocean Pole Tide.
 *
 * The file with the coefficients can be found at:
 * ftp://tai.bipm.org/iers/conv2010/chapter6/desaiscopolecoef.txt, or
 * https://iers-conventions.obspm.fr/content/chapter6/additional_info/desaiscopolecoef.txt.gz
 *
 * Note that the passed-in A_real, A_imag, B_real and B_imag matrices will be 
 * resized (at output) to hold the number of degree/order coefficients 
 * specified. Hence, at input, they can have any size.
 *
 * @param[in] fn File holding the Desai 2002 coefficients as published by IERS
 * @param[in] max_degree Max degree of coefficients to parse; should be 
 *                <= MAX_DEGREE_DESAI_2002
 * @param[in] max_order Max order of coefficients to parse; should be 
 *                <= MAX_ORDER_DESAI_2002 and <= max_degree
 * @param[out] A_real The A^{R}_{nm} coefficients of the model; will be 
 *                resized to hold the requested number of coeffs.
 * @param[out] A_imag The A^{I}_{nm} coefficients of the model; will be 
 *                resized to hold the requested number of coeffs.
 * @param[out] B_real The B^{R}_{nm} coefficients of the model; will be 
 *                resized to hold the requested number of coeffs.
 * @param[out] B_imag The B^{I}_{nm} coefficients of the model; will be 
 *                resized to hold the requested number of coeffs.
 * @return Anything other than 0, denotes an error.
 */
int parse_desai02_coeffs(
    const char *fn, int max_degree, int max_order,
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &A_real,
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &A_imag,
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &B_real,
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &B_imag) noexcept;

/** Overload of parse_desai02_coeffs with max_degree=MAX_DEGREE_DESAI_2002 
 *  and max_order=MAX_ORDER_DESAI_2002. I.e., this version will read all 
 *  available coeffs from the input file.
 *  @see parse_desai02_coeffs
 */
inline int parse_desai02_coeffs(
    const char *fn,
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &A_real,
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &A_imag,
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &B_real,
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &B_imag) noexcept {
  return parse_desai02_coeffs(fn, MAX_DEGREE_DESAI_2002, MAX_ORDER_DESAI_2002,
                              A_real, A_imag, B_real, B_imag);
}
} /* namespace pole_tide_details */

class PoleTide {
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
   *
   * TODO i do not know what time-scale t should be in
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

  /**
   */
  static CartesianCrd deformation(const MjdEpoch &t, double xp, double yp,
                                  const SphericalCrdConstView &rsta) noexcept;
}; /* class poleTide */

class OceanPoleTide {

  /* max degree of allowed expansion of spherical harmonics */
  int mmaxdegree{2};
  /* max order of allowed expansion of spherical harmonics */
  int mmaxorder{2};
  /* Stokes Coefficients; to be computed on request */
  StokesCoeffs mcs{2};
  /* Coefficients for SH expansion, according to the model of Desai 2002. See
   * also IERS 2010, Ch. 7.1.5 "Ocean pole tide loading".
   * These are used when computing the Stokes coefficients for the 
   * geopotential SH expansion, using the Desai model if and only if the 
   * mmaxdegree and mmaxorder are > 2.
   * They must be read off from the file:
   * ftp://tai.bipm.org/iers/conv2010/chapter6/desaiscopolecoef.txt
   */
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> A_real;
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> A_imag;
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> B_real;
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> B_imag;

  int parse_desai02_coeffs(const char *fn, int maxdegree,
                           int maxorder) noexcept {
    return pole_tide_details::parse_desai02_coeffs(
        fn, maxdegree, maxorder, A_real, A_imag, B_real, B_imag);
  }

public:
  int max_degree() const noexcept { return mmaxdegree; }
  int max_order() const noexcept { return mmaxorder; }
  const StokesCoeffs &stokes_coeffs() const noexcept { return mcs; }
  StokesCoeffs &stokes_coeffs() noexcept { return mcs; }
  
  OceanPoleTide(int max_degree, int max_order, const char *fn)
      : mmaxdegree(max_degree), mmaxorder(max_order),
        mcs(max_degree, max_order), A_real(max_degree, max_degree),
        A_imag(max_degree, max_degree), B_real(max_degree, max_degree),
        B_imag(max_degree, max_degree) {
    if (this->parse_desai02_coeffs(fn, max_degree, max_order)) {
      throw std::runtime_error(
          "[ERROR] Failed constructing OceanPoleTide instance\n");
    }
  }

  OceanPoleTide(int max_degree, const char *fn)
      : mmaxdegree(max_degree), mmaxorder(max_degree), mcs(max_degree),
        A_real(max_degree, max_degree), A_imag(max_degree, max_degree),
        B_real(max_degree, max_degree), B_imag(max_degree, max_degree) {
    if (this->parse_desai02_coeffs(fn, max_degree, max_degree)) {
      throw std::runtime_error(
          "[ERROR] Failed constructing OceanPoleTide instance\n");
    }
  }

  OceanPoleTide(const char *fn)
      : mmaxdegree(pole_tide_details::MAX_DEGREE_DESAI_2002),
        mmaxorder(pole_tide_details::MAX_ORDER_DESAI_2002),
        mcs(pole_tide_details::MAX_DEGREE_DESAI_2002,
            pole_tide_details::MAX_ORDER_DESAI_2002),
        A_real(pole_tide_details::MAX_DEGREE_DESAI_2002,
               pole_tide_details::MAX_DEGREE_DESAI_2002),
        A_imag(pole_tide_details::MAX_DEGREE_DESAI_2002,
               pole_tide_details::MAX_DEGREE_DESAI_2002),
        B_real(pole_tide_details::MAX_DEGREE_DESAI_2002,
               pole_tide_details::MAX_DEGREE_DESAI_2002),
        B_imag(pole_tide_details::MAX_DEGREE_DESAI_2002,
               pole_tide_details::MAX_DEGREE_DESAI_2002) {
    if (this->parse_desai02_coeffs(fn, pole_tide_details::MAX_DEGREE_DESAI_2002,
                                   pole_tide_details::MAX_ORDER_DESAI_2002)) {
      throw std::runtime_error(
          "[ERROR] Failed constructing OceanPoleTide instance\n");
    }
  }

  OceanPoleTide()
      : mmaxdegree(2), mmaxorder(2), mcs(2), A_real(0, 0), A_imag(0, 0),
        B_real(0, 0), B_imag(0, 0) {};

  /** @brief Geopotential coefficients corrections (ΔC_12 and ΔS_12) due to 
   *         ocean pole tide.
   *
   * Model implemented here is described in IERS 2010, Ch. 6.5 "Ocean 
   * pole tide". This function will only compute corrections to C21 and S21 
   * coefficients; according to IERS 2010, these account for approximately 90% 
   * of the ocean pole tide potential.
   *
   * Note that this computation (and the resulting values) implicitilly assume 
   * constants defined in IERS 2010 (e.g. for Re, GM, G, ge, etc).
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
    dS21 = -1.7232e-10 * (m2 - 0.03365 * m1);
    return 0;
  }

  /** @brief Normalized geopotential coefficients corrections Cnm, Snm due to
   *         ocean pole tide.
   *
   * Model implemented here is described in IERS 2010, Ch. 6.5 "Ocean 
   * pole tide".
   * The model implemented here, i.e. Desai 2002, assumes that the instance's 
   * matrices A_real, Aimag and B_real, B_imag are already filled with the 
   * model values up an appropriate degree and order.
   * According to IERS 2010:
   * > Approximately 90% of the variance of the ocean pole tide potential is 
   * > provided by the degree n = 2 spherical harmonic components, with the 
   * > next largest contributions provided by the degree n = 1 and n = 3 
   * > components, respectively. Expansion to spherical harmonic degree n = 10 
   * > provides approximately 99% of the variance. However, adequate 
   * > representation of the continental boundaries will require a spherical 
   * > harmonic expansion to high degree and order.
   *
   * @param[in] t Time of computation
   * @param[in] xp Pole X-coordinates at t in [arcsec]
   * @param[in] yp Pole Y-coordinates at t in [arcsec]
   * @param[in] max_degree Max degree (n) of coefficients to be computed. Note 
   *            that it should be > 0 and <= mmaxdegree
   * @param[in] max_degree Max degree (n) of coefficients to be computed. Note 
   *            that it should be > 0 and <= mmaxdegree
   * @param[in] omega Nominal mean Earth's angular velocity [rad/sec]; default 
   *            value taken from IERS 2010.
   * @param[in] G Constant of gravitation in [m^3/kg/s^2]; default value taken 
   *            from IERS 2010
   * @param[in] ge Mean equatorial gravity in [m/sec^2]; default value taken 
   *            from IERS 2010
   * @return Always zero. At output, then instance's mcs member variable will 
   * be filled with the computed coefficients, up to the the maximum 
   * degree/order specified.
   */
  int stokes_coeffs(const dso::MjdEpoch &t, double xp, double yp,
                    int max_degree, int max_order,
                    double OmegaEarth = ::iers2010::OmegaEarth,
                    double G = ::iers2010::G,
                    double ge = ::iers2010::ge) noexcept;
}; /* class oceanPoleTide */


}/* namespace dso */

#endif
