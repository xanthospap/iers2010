#ifndef __DOODSON_NUMBER_DEFINES_HPP__
#define __DOODSON_NUMBER_DEFINES_HPP__

#include <cstring>
#include <datetime/dtfund.hpp>
#include <initializer_list>
#include "geodesy/units.hpp"
#include "datetime/dtcalendar.hpp"
#include "iersc.hpp"
#include <cassert>

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

  constexpr DoodsonNumber(std::initializer_list<int> l) noexcept {
    assert(l.size() == 6);
    iar[0] = *(l.begin()+0);
    iar[1] = *(l.begin()+1);
    iar[2] = *(l.begin()+2);
    iar[3] = *(l.begin()+3);
    iar[4] = *(l.begin()+4);
    iar[5] = *(l.begin()+5);
  }

  /// @brief Transform to a Doodson-number string
  char *str(char *buf, bool use_5s_convention = false) const noexcept;

  /// @brief Get Doodson multipler at position i
  int operator()(int i) const noexcept {
#ifdef DEBUG
    assert(i>=0 && i<6);
#endif
    return iar[i];
  }

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

  /// @brief Frequency of constituent in [rad/day]
  double frequency(const double *const doodson_freq_vars) const noexcept {
    return (doodson_freq_vars[0] * iar[0] + doodson_freq_vars[1] * iar[1] +
            doodson_freq_vars[2] * iar[2] + doodson_freq_vars[3] * iar[3] +
            doodson_freq_vars[4] * iar[4] + doodson_freq_vars[5] * iar[5]) /
           dso::days_in_julian_cent;
  }

};// DoodsonNumber

/// @brief Greenwich Mean Sideral Time in [rad], range [0-2π)
/// @warning Note that this is not consistent with IAU2006[A] and IERS2010, 
/// and does not take into account ERA angle and TT time
/// TODO:
/// @see https://github.com/groops-devs/groops/blob/main/source/base/planets.cpp
/// This is not the angle iers2010::gmst??, it does not take into account
/// TT time! However, it seems that this angle is used to compute Doodson
/// arguments
/// See also the 00README_simulation.txt in COST-G benchmark
double gmst_utc(const dso::TwoPartDate &utc) noexcept;

/// @brief Fundamental (Delaunay) arguments to Doodson variables.
/// All angles are in [rad] in the range [0,2π)
///
/// @warning I THINK(!!??!!) that GMST angle here should be computed using the
///          gmst_utc function and NOT the gmst06!!
///
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
///   * [3] p  : Longitude of Moon's mean perigee [rad]
///   * [4] N' : Negative longitude of Moon's mean node [rad]
///   * [5] pl : Longitude of Sun's mean perigee [rad]
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

// TODO we need the derivative of gmst06 here!!
// Compare with what we get from TDFRPH
// The derivatives of gmst and era should be defined (once validated) in the
// iers2010 namespace
// Return values are in [rad/century]
// Input values are in [rad/century]
[[deprecated]]
inline int fundarg_derivs2doodson_derivs(const double *const fundarg_derivs,
                                         double dgmst,
                                         double *doodson) noexcept {
  const double* const f = fundarg_derivs;
  doodson[1] = f[2] + f[4];
  doodson[2] = f[2] + f[4] - f[3];
  doodson[3] = f[2] + f[4] - f[0];
  doodson[4] = -f[4];
  doodson[5] = f[2] + f[4] - f[3] - f[1];
  doodson[0] = dgmst - (f[2] + f[4]);
  return 0;
}

inline int doodson_frequency_args(const dso::TwoPartDate &tt, double *f) noexcept {
  const double t = tt.jcenturies_sinceJ2000() / 1e0;
  f[0] = dso::deg2rad(
      127037328.88553056e0 +
      (2 * 0.17696111e0 + (3 * -0.00183140e0 + 4 * 0.00008824e0 * t) * t) * t);
  f[1] = dso::deg2rad(
      4812678.81195750e0 +
      (2 * -0.14663889e0 + (3 * 0.00185140e0 + 4 * -0.00015355e0 * t) * t) * t);
  f[2] = dso::deg2rad(
      360007.69748806e0 +
      (2 * 0.03032222e0 + (3 * 0.00002000e0 + 4 * -0.00006532e0 * t) * t) * t);
  f[3] = dso::deg2rad(
      40690.13635250e0 +
      (2 * -1.03217222e0 + (3 * -0.01249168e0 + 4 * 0.00052655e0 * t) * t) * t);
  f[4] = dso::deg2rad(
      19341.36261972e0 +
      (2 * -0.20756111e0 + (3 * -0.00213942e0 + 4 * 0.00016501e0 * t) * t) * t);
  f[5] = dso::deg2rad(
      17.19457667e0 +
      (2 * 0.04568889e0 + (3 * -0.00001776e0 + 4 * -0.00003323e0 * t) * t) * t);
  for (int i=0; i<6;i ++) f[i] /= 10e0;
  return 0;
}

}//namespace dso

#endif
