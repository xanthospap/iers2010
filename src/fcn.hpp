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
 * @return Anything other than 0, denotes an error.
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
 * @return Anything other than 0, denotes an error.
 */
inline int parse_lambert_coefficients(
    const char *fn, std::vector<fcn::LambertCoefficients> &lvec) noexcept {
  return parse_lambert_coefficients(fn, TwoPartDate::min(), TwoPartDate::max(),
                                    lvec);
}
} /* namespace dso */

#endif
