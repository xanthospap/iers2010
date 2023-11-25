/** @file
 * Handle models and relevant information for modeling Free Core Nutation
 */

#ifndef __DSO_FREE_CORE_NUTATION_HPP__
#define __DSO_FREE_CORE_NUTATION_HPP__

#include "datetime/tpdate.hpp"
#include <vector>

namespace dso {
namespace fcn {
/** class to hold coefficients for Lambert's FCN model. See 
 * http://ivsopar.obspm.fr/fcn/notice.pdf
 */
struct LambertCoefficients {
  TwoPartDate t;
  double real;  /* in [microas] */
  double imag;  /* in [microas] */
  double sigma; /* in [microas] */
  LambertCoefficients(const dso::TwoPartDate &t_, double r, double i,
                      double s) noexcept
      : t(t_), real(r), imag(i), sigma(s) {}
  /* real part of coefficient, in [microas] */
  double xc() const noexcept {return real;}
  double ys() const noexcept {return -xc();}
  /* imaginary part of the coefficient, in [microas] */
  double xs() const noexcept {return imag;}
  double yc() const noexcept {return xs();}
  /* uncertainty the coefficient, in [microas] */
  double sx() const noexcept {return sigma;}
  double sy() const noexcept {return sx();}
}; /* LambertCoefficients */

/** A struct to hold FCN results, i.e. (X,Y) of CIP in GCRS */
struct FcnResult {
  double xcip; /* in [microas] */
  double ycip; /* in [microas] */
  double sxcip; /* in [microas] */
  double sycip; /* in [microas] */
  FcnResult(double xc, double yc, double sx, double sy) noexcept
      : xcip(xc), ycip(yc), sxcip(sx), sycip(sy) {}
}; /* FcnResult */

std::vector<LambertCoefficients> load_iers10_table52c() noexcept;
} /* namespace fcn */

/** @brief Parse a Lambert coefficients (ascii) file for FCN modeling.
 *
 * These files are distributed at http://ivsopar.obspm.fr/fcn/ to model FCN 
 * using Lambert's model, see http://ivsopar.obspm.fr/fcn/notice.pdf.
 * 
 * @param[in] fn Name of the coefficients file
 * @param[in] from Parse coefficients using this starting date (inclusive)
 * @param[in] to Parse coefficients using this ending date (exclusive)
 * @param[in] lvec At output, a vector of LambertCoefficients parsed from 
 *            the input file, matching the interval [from, to)
 * @return Anything other than zero denotes an error.
 * epoch.
 */
int parse_lambert_coefficients(
    const char *fn, const TwoPartDate &from, const TwoPartDate &to,
    std::vector<fcn::LambertCoefficients> &lvec) noexcept;

/** @brief Parse a Lambert coefficients (ascii) file for FCN modeling.
 *
 * These files are distributed at http://ivsopar.obspm.fr/fcn/ to model FCN 
 * using Lambert's model, see http://ivsopar.obspm.fr/fcn/notice.pdf.
 * 
 * @param[in] fn Name of the coefficients file
 * @param[in] lvec At output, a vector of LambertCoefficients parsed from 
 *            the input file.
 * @return Anython other than zero denotes an error.
 */
inline int parse_lambert_coefficients(
    const char *fn, std::vector<fcn::LambertCoefficients> &lvec) noexcept {
      const auto from = TwoPartDate::min();
      const auto to = TwoPartDate::max();
  return parse_lambert_coefficients(fn, from, to,
                                    lvec);
}

/** @brief Compute CIP offset (and uncertainties) in GCRS due to FCN.
 *
 * This is an empirical model due to Lambert (see [2]) to model FCN. It 
 * computes quantities to be added to the X, Y coordinates of the CIP in the 
 * GCRS to account for the FCN effect (see [1]).
 *
 * Note that this function is NOT the one distributed by the IERS. It uses 
 * a user-defined data set (\p lvec), where the IERS-distributed FCNNUT.f 
 * routine has a hard-coded constant table for the coefficients of the 
 * empirical model (of the retrograde FCN) during the interval 1984-2012 
 * (see the included file Table5.2.c).
 *
 * References:
 * [1] IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). (IERS 
 * Technical Note ; 36) Frankfurt am Main: Verlag des Bundesamts für 
 * Kartographie und Geodäsie, 2010. 179 pp., ISBN 3-89888-989-6, 
 * Section 5.5.5
 * [2] Lambert, S., 2007, “Empirical modeling of the retrograde Free Core 
 * Nutation,” available at ftp://hpiers.obspm.fr/eop-pc/models/fcn/notice.pdf.
 * [3] http://ivsopar.obspm.fr/fcn/index.html
 */
fcn::FcnResult
lambert_fcn(const TwoPartDate &t,
            const std::vector<fcn::LambertCoefficients> &lvec) noexcept;

/** @brief Compute CIP offset (and uncertainties) in GCRS due to FCN.
 *
 * This is an empirical model due to Lambert (see [2]) to model FCN. It 
 * computes quantities to be added to the X, Y coordinates of the CIP in the 
 * GCRS to account for the FCN effect (see [1]).
 *
 * This function is an alternate implementation of the IERS-distributed 
 * FCNNUT.f routine. It uses the data provided by IERS based on the Lambert 
 * 2007 model (see the included file Table5.2.c).
 *
 * References:
 * [1] IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). (IERS 
 * Technical Note ; 36) Frankfurt am Main: Verlag des Bundesamts für 
 * Kartographie und Geodäsie, 2010. 179 pp., ISBN 3-89888-989-6, 
 * Section 5.5.5
 * [2] Lambert, S., 2007, “Empirical modeling of the retrograde Free Core 
 * Nutation,” available at ftp://hpiers.obspm.fr/eop-pc/models/fcn/notice.pdf.
 */
inline fcn::FcnResult lambert_fcn(const TwoPartDate &t) noexcept {
  return lambert_fcn(t, fcn::load_iers10_table52c());
}
} /* namespace dso */

#endif
