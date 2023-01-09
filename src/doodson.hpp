#ifndef __DOODSON_NUMBER_DEFINES_HPP__
#define __DOODSON_NUMBER_DEFINES_HPP__

#include <cstring>
#include <datetime/dtfund.hpp>
#include "geodesy/units.hpp"
#include "ggdatetime/dtcalendar.hpp"

namespace dso {

/// @class Doodson
/// @ref "Indexing and argument conventions for tides", Richard Ray (GSFC), 
///       2017, available at: 
/// https://ivscc.gsfc.nasa.gov/hfeop_wg/memos/memo-conventions_Ray_2017Dec10.pdf
class DoodsonNumber {
private:
  /// the six first arguemnts are the Doodson multipliers for Doodson
  /// arguments, [τ, s, h, p, N', p_s]
  /// the seventh term, (iar[6]) is an integer multiple of 90[deg] (see
  /// Ray, 2017)
  /// @warning We store here the "multipliers" of the Doodson variables,
  ///          without the +/-5 convention oftenly used (following Ray, 2017)
  int iar[7] = {0};

public:
  explicit DoodsonNumber(const int *ar = nullptr) {
    if (ar)
      std::memcpy(iar, ar, sizeof(int) * 6);
  }

  /// @brief Transform to a Doodson-number string
  char *str(char *buf, bool use_5s_convention = false) const noexcept;

  /// @brief Equality comparisson
  bool operator==(const DoodsonNumber &other) const noexcept {
    return (iar[0] == other.iar[0] && iar[1] == other.iar[1] &&
            iar[2] == other.iar[2] && iar[3] == other.iar[3] &&
            iar[4] == other.iar[4] && iar[5] == other.iar[5] &&
            iar[6] == other.iar[6]);
  }

  /// @brief Inequality comparisson
  bool operator!=(const DoodsonNumber &other) const noexcept {
    return !(this->operator==(other));
  }

  /// @brief Compute θ_f = Σ(β_i * n_i), i=1,..,6
  double phase(const double *const doodson_vars) const noexcept {
    return dso::anp(doodson_vars[0] * iar[0] + doodson_vars[1] * iar[1] +
                    doodson_vars[2] * iar[2] + doodson_vars[3] * iar[3] +
                    doodson_vars[4] * iar[4] + doodson_vars[5] * iar[5]);
  }
  double phase(const dso::TwoPartDate &tt_mjd,
                     const dso::TwoPartDate &ut1_mjd) const noexcept;

  static double *doodson_freq_vars(const dso::TwoPartDate &tt_mjd,
                                   double *args) noexcept;
  double frequency(const double *const doodson_freq_vars) const noexcept {
    return (doodson_freq_vars[0] * iar[0] + doodson_freq_vars[1] * iar[1] +
            doodson_freq_vars[2] * iar[2] + doodson_freq_vars[3] * iar[3] +
            doodson_freq_vars[4] * iar[4] + doodson_freq_vars[5] * iar[5]) /
           dso::days_in_julian_cent;
  }
  double frequency(const dso::TwoPartDate &tt_mjd) noexcept {
    double fargs[6];
    DoodsonNumber::doodson_freq_vars(tt_mjd, fargs);
    return frequency(fargs);
  }
};// DoodsonNumber

/// @brief Fundamental (Delaunay) arguments to Doodson variables.
/// All angles are in [rad] in the range [0,2π)
/// @param[in] fundarg Fundamental (Delaunay) arguments, in the order
///             [l, lp, f, d, Ω], see notes.
/// @param[out] doodson Corresponding Doodson variables, in the order
///             [τ, s, h, p, N', ps]
/// @note Explanation of symbols used:
///   * [0] l  : Mean anomaly of the Moon [rad]
///   * [1] lp : Mean anomaly of the Sun [rad]
///   * [2] f  : L - Ω [rad]
///   * [3] d  : Mean elongation of the Moon from the Sun [rad]
///   * [4] Ω  : Mean longitude of the ascending node of the Moon [rad]
///
///   * [0] τ  : GMST + π - s
///   * [1] s  : Moon's mean longitude [rad]
///   * [2] h  : Sun's mean longitude [rad]
///   * [3] p  : Longitude of Moon's mean perigee
///   * [4] N' : Negative longitude of Moon's mean node
///   * [5] pl : Longitude of Sun's mean perigee
inline int fundarg2doodson(const double *const fundarg, double gmst,
                           double *doodson) noexcept {
  doodson[1] = dso::anp(fundarg[2] + fundarg[4]);
  doodson[2] = dso::anp(fundarg[2] + fundarg[4] - fundarg[3]);
  doodson[3] = dso::anp(fundarg[2] + fundarg[4] - fundarg[0]);
  doodson[4] = dso::anp(-fundarg[4]);
  doodson[5] = dso::anp(fundarg[2] + fundarg[4] - fundarg[3] - fundarg[1]);
  doodson[0] = dso::anp(gmst + dso::DPI - doodson[1]);
  return 0;
}
}//namespace dso

#endif
