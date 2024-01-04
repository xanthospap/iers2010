/** @file
 * Define a DoodsonConstituent class to represent tidal constituents via
 * Doodson numbers.
 * References:
 * [1]
 * https://ivscc.gsfc.nasa.gov/hfeop_wg/memos/memo-conventions_Ray_2017Dec10.pdf
 */

#ifndef __DOODSON_NUMBER_DEFINES_HPP__
#define __DOODSON_NUMBER_DEFINES_HPP__

#include "datetime/calendar.hpp"
#include "geodesy/units.hpp"
#include <array>
#include <cassert>
#include <cstring>
#include <initializer_list>
#include <stdexcept>

namespace dso {

/* @class
 * @ref "Indexing and argument conventions for tides", Richard Ray (GSFC),
 *        2017, available at [1]
 */
class DoodsonConstituent {
private:
  /* The six first arguemnts are the Doodson multipliers for Doodson
   * arguments, [τ, s, h, p, N', p_s]
   * The seventh term, (pifac) is an integer multiple of 90[deg] (see
   * Ray, 2017)
   *
   * @warning We store here the "multipliers" of the Doodson variables,
   * without the +/-5 convention oftenly used (following Ray, 2017)
   */
  int iar[6] = {0};
  double pifac{0e0};

public:
  /* @brief Constructor given (optionally) an int array. Default values (if
   *        no array is given) are 0
   * @param[in] ar An array of 6 integers, interpreted as:
   *  ar[0] multiplier for τ   
   *  ar[1] multiplier for s   
   *  ar[2] multiplier for h   
   *  ar[3] multiplier for p   
   *  ar[4] multiplier for N'  
   *  ar[5] multiplier for p_s 
   */
  explicit DoodsonConstituent(const int *ar = nullptr) {
    if (ar)
      std::memcpy(iar, ar, sizeof(int) * 6);
  }

  /* @brief Constructor from initializer list (of ints)
   *        i.e. enables: DoodsonConstituent d(1,2,3,4,5,6);
   * @param[in] l An initializer list of 6 ints, interpreted as:
   *   l[0] multiplier for τ   
   *   l[1] multiplier for s   
   *   l[2] multiplier for h   
   *   l[3] multiplier for p   
   *   l[4] multiplier for N'  
   *   l[5] multiplier for p_s 
   */
  constexpr DoodsonConstituent(std::initializer_list<int> l) noexcept {
    assert(l.size() == 6);
    iar[0] = *(l.begin() + 0);
    iar[1] = *(l.begin() + 1);
    iar[2] = *(l.begin() + 2);
    iar[3] = *(l.begin() + 3);
    iar[4] = *(l.begin() + 4);
    iar[5] = *(l.begin() + 5);
  }

  /* @brief Constructor from initializer list (of ints) and π factor
   *        i.e. enables: DoodsonConstituent d(1,2,3,4,5,6,0.5);
   * @param[in] a0 : τ
   * @param[in] a1 : s
   * @param[in] a2 : h
   * @param[in] a3 : p
   * @param[in] a4 : N'
   * @param[in] a5 : p_s
   * @param[in] pifactor : π-factor
   */
  constexpr DoodsonConstituent(int a0, int a1, int a2, int a3, int a4, int a5,
                               double pifactor) noexcept {
    iar[0] = a0;
    iar[1] = a1;
    iar[2] = a2;
    iar[3] = a3;
    iar[4] = a4;
    iar[5] = a5;
    pifac = pifactor;
  }

  /* @brief Transform to a Doodson-number string as: "%1d%2d%2d.%2d%2d%2d"
   * @param[out] buf A char buffer of size at least 13 where the instance
   *            will be written to.
   * @param[in] use_5s_convention Add +5 to all integers after the first
   *            place (i.e. write the string as: τ,s+5,h+5,p+5,N'+5,ps+5)
   * @return A pointer to the start of the string, i.e. to buf
   */
  char *str(char *buf, bool use_5s_convention = false) const noexcept;

  /* @brief Get Doodson multipler/integer at position i */
  int operator()(int i) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 6);
#endif
    return iar[i];
  }

  /* @brief Get Doodson multipler at position i */
  int &operator()(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 6);
#endif
    return iar[i];
  }

  /* @brief Equality comparisson */
  bool operator==(const DoodsonConstituent &other) const noexcept {
    return (iar[0] == other.iar[0] && iar[1] == other.iar[1] &&
            iar[2] == other.iar[2] && iar[3] == other.iar[3] &&
            iar[4] == other.iar[4] && iar[5] == other.iar[5] &&
            pifac == other.pifac);
  }

  /* @brief Inequality comparisson */
  bool operator!=(const DoodsonConstituent &other) const noexcept {
    return !(this->operator==(other));
  }

  /* @brief Check if two tidal constituents are of the same "spicies", i.e.
   *        share k1 values (see [1]).
   */
  constexpr bool same_species(const DoodsonConstituent &other) const noexcept {
    return iar[0] == other.iar[0];
  }

  /* @brief Check if two tidal constituents are of the same "group", i.e.
   *        share k1 and k2 values (see [1]).
   */
  constexpr bool same_group(const DoodsonConstituent &other) const noexcept {
    return (this->same_species(other) && (iar[1] == other.iar[1]));
  }

  /* @brief Check if two tidal constituents are of the same "constituent", i.e.
   *        share k1, k2 and k3 values (see [1]).
   */
  constexpr bool
  same_constituent(const DoodsonConstituent &other) const noexcept {
    return (this->same_group(other) && (iar[2] == other.iar[2]));
  }

  /** @brief Compute the 'angular' argument (ARGUMENT) 
   *
   * The argument is often designated as ARGUMENT in the iERS standards and 
   * often is also called 'phase'. 
   * The computation formula is: 
   *           θ_f = Σ(β_i * n_i), i=1,..,6
   * where:
   * β_i are the [τ, s, h, p, N', pl] arguments (e.g. derived from Delaunay 
   * arguments plus GMST) and n_i are the Doodson multipliers.
   *
   * @return Angular argument in [rad], within the range [0,2π).
   */
  double
  argument(const double *const doodson_arguments) const noexcept {
    const double *__restrict__ f = doodson_arguments;
    return dso::anp(f[0] * iar[0] + f[1] * iar[1] + f[2] * iar[2] +
                    f[3] * iar[3] + f[4] * iar[4] + f[5] * iar[5]);
  }

}; /* DoodsonConstituent */

/* @brief Fundamental (Delaunay) arguments to Doodson variables.
 * All angles are in [rad] in the range [0,2π)
 *
 * @param[in] fundarg Fundamental (Delaunay) arguments, in the order
 *             [l, lp, f, d, Ω], see notes.
 * @param[out] doodson Corresponding Doodson variables, in the order
 *             [τ, s, h, p, N', ps]
 * 
 * @note Explanation of symbols used:
 *   * [0] l  : Mean anomaly of the Moon [rad]
 *   * [1] lp : Mean anomaly of the Sun [rad]
 *   * [2] f  : L - Ω [rad]
 *   * [3] d  : Mean elongation of the Moon from the Sun [rad]
 *   * [4] Ω  : Mean longitude of the ascending node of the Moon [rad]
 *
 *   * [0] τ  : GMST + π - s
 *   * [1] s  : Moon's mean longitude [rad]
 *   * [2] h  : Sun's mean longitude [rad]
 *   * [3] p  : Longitude of Moon's mean perigee [rad]
 *   * [4] N' : Negative longitude of Moon's mean node [rad]
 *   * [5] pl : Longitude of Sun's mean perigee [rad]
 */
inline double *delaunay2doodson(const double *const fundarg, double gmst,
                           double *doodson) noexcept {
  doodson[1] = dso::anp(fundarg[2] + fundarg[4]);
  doodson[2] = dso::anp(fundarg[2] + fundarg[4] - fundarg[3]);
  doodson[3] = dso::anp(fundarg[2] + fundarg[4] - fundarg[0]);
  doodson[4] = dso::anp(-fundarg[4]);
  doodson[5] = dso::anp(fundarg[2] + fundarg[4] - fundarg[3] - fundarg[1]);
  doodson[0] = dso::anp(gmst + dso::DPI - doodson[1]);
  return doodson;
}

namespace detail {
struct TidalConstituentsArrayEntry {
  DoodsonConstituent _d;
  double _per; /* period in [deg/hour] */
  const char _n[8]; /* name */
}; /* TidalConstituentsArrayEntry */
} /* namespace detail */

/** Given a tidal constituent name, return its details */
detail::TidalConstituentsArrayEntry get_wave(const char *name);

constexpr static std::array<detail::TidalConstituentsArrayEntry, 166>
    TidalConstituentsArray = {{
        {/*055565*/ DoodsonConstituent(0, +0, +0, +0, +1, +0, +1.0e+00),
         0.00221000, ""},
        {/*055575*/ DoodsonConstituent(0, +0, +0, +0, +2, +0, +0.0e+00),
         0.00441000, ""},
        {/*056554*/ DoodsonConstituent(0, +0, +1, +0, +0, -1, +0.0e+00),
         0.04107000, ""},
        {/*056555*/ DoodsonConstituent(0, +0, +1, +0, +0, +0, +0.0e+00),
         0.04106864, ""},
        {/*057555*/ DoodsonConstituent(0, +0, +2, +0, +0, +0, +0.0e+00),
         0.08213728, ""},
        {/*057565*/ DoodsonConstituent(0, +0, +2, +0, +1, +0, -9.0e+00),
         0.08434000, ""}, /* TODO unknown π multiplier  */
        {/*058554*/ DoodsonConstituent(0, +0, +3, +0, +0, -1, +0.0e+00),
         0.12320000, ""},
        {/*063655*/ DoodsonConstituent(0, +1, -2, +1, +0, +0, +0.0e+00),
         0.47152109, ""},
        {/*065445*/ DoodsonConstituent(0, +1, +0, -1, -1, +0, -9.0e+00),
         0.54217000, ""}, /* TODO unknown π multiplier  */
        {/*065455*/ DoodsonConstituent(0, +1, +0, -1, +0, +0, +0.0e+00),
         0.54437470, ""},
        {/*065465*/ DoodsonConstituent(0, +1, +0, -1, +1, +0, -9.0e+00),
         0.54658000, ""}, /* TODO unknown π multiplier  */
        {/*065655*/ DoodsonConstituent(0, +1, +0, +1, +0, +0, -9.0e+00),
         0.55366000, ""}, /* TODO unknown π multiplier  */
        {/*073555*/ DoodsonConstituent(0, +2, -2, +0, +0, +0, +0.0e+00),
         1.01589579, ""},
        {/*075355*/ DoodsonConstituent(0, +2, +0, -2, +0, +0, -9.0e+00),
         1.08875000, ""}, /* TODO unknown π multiplier  */
        {/*075555*/ DoodsonConstituent(0, +2, +0, +0, +0, +0, +0.0e+00),
         1.09803306, ""},
        {/*075565*/ DoodsonConstituent(0, +2, +0, +0, +1, +0, -9.0e+00),
         1.10024000, ""}, /* TODO unknown π multiplier  */
        {/*075575*/ DoodsonConstituent(0, +2, +0, +0, +2, +0, -9.0e+00),
         1.10245000, ""}, /* TODO unknown π multiplier  */
        {/*083655*/ DoodsonConstituent(0, +3, -2, +1, +0, +0, +0.0e+00),
         1.56955415, ""},
        {/*085455*/ DoodsonConstituent(0, +3, +0, -1, +0, +0, +0.0e+00),
         1.64240776, ""},
        {/*085465*/ DoodsonConstituent(0, +3, +0, -1, +1, +0, -9.0e+00),
         1.64462000, ""}, /* TODO unknown π multiplier  */
        {/*093555*/ DoodsonConstituent(0, +4, -2, +0, +0, +0, +0.0e+00),
         2.11394000, ""},
        {/*095355*/ DoodsonConstituent(0, +4, +0, -2, +0, +0, -0.0e+00),
         2.18678246, ""},
        {/*105555*/ DoodsonConstituent(1, -5, +0, +0, +0, +0, +2.5e+00),
         11.74696945, ""},
        {/*125755*/ DoodsonConstituent(1, -3, +0, +2, +0, +0, +5.0e-01),
         12.85428618, ""},
        {/*127555*/ DoodsonConstituent(1, -3, +2, +0, +0, +0, +5.0e-01),
         12.92713979, ""},
        {/*135645*/ DoodsonConstituent(1, -2, +0, +1, -1, +0, -9.0e+00),
         13.39645000, ""}, /* TODO unknown π multiplier  */
        {/*135655*/ DoodsonConstituent(1, -2, +0, +1, +0, +0, +5.0e-01),
         13.39866088, ""},
        {/*137455*/ DoodsonConstituent(1, -2, +2, -1, +0, +0, +5.0e-01),
         13.47151449, ""},
        {/*145545*/ DoodsonConstituent(1, -1, +0, +0, -1, +0, -9.0e+00),
         13.94083000, ""}, /* TODO unknown π multiplier  */
        {/*145555*/ DoodsonConstituent(1, -1, +0, +0, +0, +0, +5.0e-01),
         13.94303557, ""},
        {/*146555*/ DoodsonConstituent(1, -1, +1, +0, +0, +0, -0.0e+00),
         13.98410421, ""},
        {/*147555*/ DoodsonConstituent(1, -1, +2, +0, +0, +0, -5.0e-01),
         14.02517285, ""},
        {/*153655*/ DoodsonConstituent(1, +0, -2, +1, +0, +0, -9.0e+00),
         14.41456000, ""}, /* TODO unknown π multiplier  */
        {/*155445*/ DoodsonConstituent(1, +0, +0, -1, -1, +0, -9.0e+00),
         14.48520000, ""}, /* TODO unknown π multiplier  */
        {/*155455*/ DoodsonConstituent(1, +0, +0, -1, +0, +0, -9.0e+00),
         14.48741000, ""}, /* TODO unknown π multiplier  */
        {/*155655*/ DoodsonConstituent(1, +0, +0, +1, +0, +0, -5.0e-01),
         14.49669394, ""},
        {/*155665*/ DoodsonConstituent(1, +0, +0, +1, +1, +0, -9.0e+00),
         14.49890000, ""}, /* TODO unknown π multiplier  */
        {/*157455*/ DoodsonConstituent(1, +0, +2, -1, +0, +0, -5.0e-01),
         14.56954755, ""},
        {/*157465*/ DoodsonConstituent(1, +0, +2, -1, +1, +0, -9.0e+00),
         14.57176000, ""}, /* TODO unknown π multiplier  */
        {/*162556*/ DoodsonConstituent(1, +1, -3, +0, +0, +1, +5.0e-01),
         14.91786468, ""},
        {/*163545*/ DoodsonConstituent(1, +1, -2, +0, -1, +0, -9.0e+00),
         14.95673000, ""}, /* TODO unknown π multiplier  */
        {/*163555*/ DoodsonConstituent(1, +1, -2, +0, +0, +0, +5.0e-01),
         14.95893136, "P1"},
        {/*164554*/ DoodsonConstituent(1, +1, -1, +0, +0, -1, -9.0e+00),
         15.00000000, ""}, /* TODO unknown π multiplier  */
        {/*164555*/ DoodsonConstituent(1, +1, -1, +0, +0, +0, +1.0e+00),
         15.00000000, ""},
        {/*164556*/ DoodsonConstituent(1, +1, -1, +0, +0, +1, +1.0e+00),
         15.00000000, "S1"},
        {/*165345*/ DoodsonConstituent(1, +1, +0, -2, -1, +0, -9.0e+00),
         15.02958000, ""}, /* TODO unknown π multiplier  */
        {/*165535*/ DoodsonConstituent(1, +1, +0, +0, -2, +0, -9.0e+00),
         15.03665000, ""}, /* TODO unknown π multiplier  */
        {/*165545*/ DoodsonConstituent(1, +1, +0, +0, -1, +0, -9.0e+00),
         15.03886000, ""}, /* TODO unknown π multiplier  */
        {/*165555*/ DoodsonConstituent(1, +1, +0, +0, +0, +0, -5.0e-01),
         15.04106864, "K1"},
        {/*165565*/ DoodsonConstituent(1, +1, +0, +0, +1, +0, -9.0e+00),
         15.04328000, ""}, /* TODO unknown π multiplier  */
        {/*165575*/ DoodsonConstituent(1, +1, +0, +0, +2, +0, -9.0e+00),
         15.04548000, ""}, /* TODO unknown π multiplier  */
        {/*166455*/ DoodsonConstituent(1, +1, +1, -1, +0, +0, -9.0e+00),
         15.07749000, ""}, /* TODO unknown π multiplier  */
        {/*166544*/ DoodsonConstituent(1, +1, +1, +0, -1, -1, -9.0e+00),
         15.07993000, ""}, /* TODO unknown π multiplier  */
        {/*166554*/ DoodsonConstituent(1, +1, +1, +0, +0, -1, -5.0e-01),
         15.08213532, ""},
        {/*166556*/ DoodsonConstituent(1, +1, +1, +0, +0, +1, -9.0e+00),
         15.08214000, ""}, /* TODO unknown π multiplier  */
        {/*166564*/ DoodsonConstituent(1, +1, +1, +0, +1, -1, -9.0e+00),
         15.08434000, ""}, /* TODO unknown π multiplier  */
        {/*167355*/ DoodsonConstituent(1, +1, +2, -2, +0, +0, -9.0e+00),
         15.11392000, ""}, /* TODO unknown π multiplier  */
        {/*167365*/ DoodsonConstituent(1, +1, +2, -2, +1, +0, -9.0e+00),
         15.11613000, ""}, /* TODO unknown π multiplier  */
        {/*167555*/ DoodsonConstituent(1, +1, +2, +0, +0, +0, -5.0e-01),
         15.12320592, ""},
        {/*167565*/ DoodsonConstituent(1, +1, +2, +0, +1, +0, -9.0e+00),
         15.12542000, ""}, /* TODO unknown π multiplier  */
        {/*168554*/ DoodsonConstituent(1, +1, +3, +0, +0, -1, -9.0e+00),
         15.16427000, ""}, /* TODO unknown π multiplier  */
        {/*173655*/ DoodsonConstituent(1, +2, -2, +1, +0, +0, -5.0e-01),
         15.51258973, ""},
        {/*173665*/ DoodsonConstituent(1, +2, -2, +1, +1, +0, -9.0e+00),
         15.51480000, ""}, /* TODO unknown π multiplier  */
        {/*175445*/ DoodsonConstituent(1, +2, +0, -1, -1, +0, -9.0e+00),
         15.58323000, ""}, /* TODO unknown π multiplier  */
        {/*175455*/ DoodsonConstituent(1, +2, +0, -1, +0, +0, -5.0e-01),
         15.58544334, ""},
        {/*175465*/ DoodsonConstituent(1, +2, +0, -1, +1, +0, -9.0e+00),
         15.58765000, ""}, /* TODO unknown π multiplier  */
        {/*183555*/ DoodsonConstituent(1, +3, -2, +0, +0, +0, -5.0e-01),
         16.05696443, ""},
        {/*185355*/ DoodsonConstituent(1, +3, +0, -2, +0, +0, -9.0e+00),
         16.12989000, ""}, /* TODO unknown π multiplier  */
        {/*185555*/ DoodsonConstituent(1, +3, +0, +0, +0, +0, -5.0e-01),
         16.13910170, ""},
        {/*185565*/ DoodsonConstituent(1, +3, +0, +0, +1, +0, -9.0e+00),
         16.14131000, ""}, /* TODO unknown π multiplier  */
        {/*185575*/ DoodsonConstituent(1, +3, +0, +0, +2, +0, -9.0e+00),
         16.14352000, ""}, /* TODO unknown π multiplier  */
        {/*195455*/ DoodsonConstituent(1, +4, +0, -1, +0, +0, -5.0e-01),
         16.68347640, ""},
        {/*195465*/ DoodsonConstituent(1, +4, +0, -1, +1, +0, -9.0e+00),
         16.68569000, ""}, /* TODO unknown π multiplier  */
        {/*209655*/ DoodsonConstituent(2, -5, +4, +1, +0, +0, -0.0e+00),
         26.40793794, ""},
        {/*217755*/ DoodsonConstituent(2, -4, +2, +2, +0, +0, -0.0e+00),
         26.87945903, ""},
        {/*219555*/ DoodsonConstituent(2, -4, +4, +0, +0, +0, -0.0e+00),
         26.95231264, ""},
        {/*219755*/ DoodsonConstituent(2, -4, +4, +2, +0, +0, -0.0e+00),
         26.96159631, ""},
        {/*225655*/ DoodsonConstituent(2, -3, +0, +1, +0, +0, +1.0e+00),
         27.34169645, ""},
        {/*225855*/ DoodsonConstituent(2, -3, +0, +3, +0, +0, +0.0e+00),
         0.00000000, ""}, /* TODO unknown period  */
        {/*227655*/ DoodsonConstituent(2, -3, +2, +1, +0, +0, +0.0e+00),
         27.42383373, ""},
        {/*229455*/ DoodsonConstituent(2, -3, +4, -1, +0, +0, -0.0e+00),
         27.49668734, ""},
        {/*229655*/ DoodsonConstituent(2, -3, +4, +1, +0, +0, -0.0e+00),
         27.50597101, ""},
        {/*233555*/ DoodsonConstituent(2, -2, -2, +0, +0, +0, -0.0e+00),
         27.80393387, ""},
        {/*235555*/ DoodsonConstituent(2, -2, +0, +0, +0, +0, -0.0e+00),
         27.88607115, ""},
        {/*235755*/ DoodsonConstituent(2, -2, +0, +2, +0, +0, +0.0e+00),
         27.89535482, ""},
        {/*237555*/ DoodsonConstituent(2, -2, +2, +0, +0, +0, +0.0e+00),
         27.96820843, ""},
        {/*243655*/ DoodsonConstituent(2, -1, -2, +1, +0, +0, -0.0e+00),
         28.35759224, ""},
        {/*245655*/ DoodsonConstituent(2, -1, +0, +1, +0, +0, +0.0e+00),
         28.43972952, "N2"},
        {/*247455*/ DoodsonConstituent(2, -1, +2, -1, +0, +0, +0.0e+00),
         28.51258312, ""},
        {/*253555*/ DoodsonConstituent(2, +0, -2, +0, +0, +0, -0.0e+00),
         28.90196694, ""},
        {/*253755*/ DoodsonConstituent(2, +0, -2, +2, +0, +0, +1.0e+00),
         0.00000000, ""}, /* TODO unknown period  */
        {/*254555*/ DoodsonConstituent(2, +0, -1, +0, +0, +0, +5.0e-01),
         28.94303557, ""},
        {/*254556*/ DoodsonConstituent(2, +0, -1, +0, +0, +1, +1.0e+00),
         0.00000000, ""}, /* TODO unknown period  */
        {/*255555*/ DoodsonConstituent(2, +0, +0, +0, +0, +0, +0.0e+00),
         28.98410421, "M2"},
        {/*256554*/ DoodsonConstituent(2, +0, +1, +0, +0, -1, +0.0e+00),
         0.00000000, ""}, /* TODO unknown period  */
        {/*256555*/ DoodsonConstituent(2, +0, +1, +0, +0, +0, -5.0e-01),
         29.02517285, ""},
        {/*257555*/ DoodsonConstituent(2, +0, +2, +0, +0, +0, +0.0e+00),
         29.06624149, ""},
        {/*263655*/ DoodsonConstituent(2, +1, -2, +1, +0, +0, +1.0e+00),
         29.45562530, ""},
        {/*265455*/ DoodsonConstituent(2, +1, +0, -1, +0, +0, +1.0e+00),
         29.52847891, "L2"},
        {/*265655*/ DoodsonConstituent(2, +1, +0, +1, +0, +0, -0.0e+00),
         29.53776258, ""},
        {/*271557*/ DoodsonConstituent(2, +2, -4, +0, +0, +2, +0.0e+00),
         0.00000000, ""}, /* TODO unknown period  */
        {/*272556*/ DoodsonConstituent(2, +2, -3, +0, +0, +1, +0.0e+00),
         29.95893332, "T2"},
        {/*273555*/ DoodsonConstituent(2, +2, -2, +0, +0, +0, +0.0e+00),
         30.00000000, "S2"},
        {/*274554*/ DoodsonConstituent(2, +2, -1, +0, +0, -1, +1.0e+00),
         30.04106668, "R2"},
        {/*275555*/ DoodsonConstituent(2, +2, +0, +0, +0, +0, +0.0e+00),
         30.08213728, ""},
        {/*283455*/ DoodsonConstituent(2, +3, -2, -1, +0, +0, -0.0e+00),
         30.54437470, ""},
        {/*283655*/ DoodsonConstituent(2, +3, -2, +1, +0, +0, +0.0e+00),
         0.00000000, ""}, /* TODO unknown period  */
        {/*285455*/ DoodsonConstituent(2, +3, +0, -1, +0, +0, +0.0e+00),
         30.62651197, ""},
        {/*291555*/ DoodsonConstituent(2, +4, -4, +0, +0, +0, -0.0e+00),
         31.01589579, ""},
        {/*293555*/ DoodsonConstituent(2, +4, -2, +0, +0, +0, -0.0e+00),
         31.09803306, ""},
        {/*301455*/ DoodsonConstituent(3, -5, -4, -1, +0, +0, -0.0e+00),
         31.56027048, ""},
        {/*309555*/ DoodsonConstituent(3, -5, +4, +0, +0, +0, -0.0e+00),
         32.03179157, ""},
        {/*335655*/ DoodsonConstituent(3, -2, +0, +1, +0, +0, +5.0e-01),
         42.38276509, ""},
        {/*345555*/ DoodsonConstituent(3, -1, +0, +0, +0, +0, +5.0e-01),
         42.92713979, ""},
        {/*355555*/ DoodsonConstituent(3, +0, +0, +0, +0, +0, +1.0e+00),
         43.47615632, ""},
        {/*363555*/ DoodsonConstituent(3, +1, -2, +0, +0, +0, +5.0e-01),
         43.94303558, ""},
        {/*364555*/ DoodsonConstituent(3, +1, -1, +0, +0, +0, -0.0e+00),
         43.98410421, ""},
        {/*365555*/ DoodsonConstituent(3, +1, +0, +0, +0, +0, -5.0e-01),
         44.02517285, ""},
        {/*381555*/ DoodsonConstituent(3, +3, -4, +0, +0, +0, +1.0e+00),
         44.95893136, "T3"},
        {/*382555*/ DoodsonConstituent(3, +3, -3, +0, +0, +0, +1.0e+00),
         45.00000000, "S3"},
        {/*383555*/ DoodsonConstituent(3, +3, -2, +0, +0, +0, +1.0e+00),
         45.04106864, "R3"},
        {/*385555*/ DoodsonConstituent(3, +3, +0, +0, +0, +0, -1.5e+00),
         45.12320592, ""},
        {/*427655*/ DoodsonConstituent(4, -3, +2, +1, +0, +0, -0.0e+00),
         56.40793794, ""},
        {/*435755*/ DoodsonConstituent(4, -2, +0, +2, +0, +0, +0.0e+00),
         56.87945903, ""},
        {/*437555*/ DoodsonConstituent(4, -2, +2, +0, +0, +0, -0.0e+00),
         56.95231264, ""},
        {/*445655*/ DoodsonConstituent(4, -1, +0, +1, +0, +0, +0.0e+00),
         57.42383373, ""},
        {/*447455*/ DoodsonConstituent(4, -1, +2, -1, +0, +0, -0.0e+00),
         57.49668734, ""},
        {/*455555*/ DoodsonConstituent(4, +0, +0, +0, +0, +0, +0.0e+00),
         57.96820843, ""},
        {/*463655*/ DoodsonConstituent(4, +1, -2, +1, +0, +0, -0.0e+00),
         58.43972952, ""},
        {/*465455*/ DoodsonConstituent(4, +1, +0, -1, +0, +0, +1.0e+00),
         58.51258312, ""},
        {/*465655*/ DoodsonConstituent(4, +1, +0, +1, +0, +0, -0.0e+00),
         58.52186679, ""},
        {/*473555*/ DoodsonConstituent(4, +2, -2, +0, +0, +0, +0.0e+00),
         58.98410421, ""},
        {/*475555*/ DoodsonConstituent(4, +2, +0, +0, +0, +0, -0.0e+00),
         59.06624149, ""},
        {/*483455*/ DoodsonConstituent(4, +3, -2, -1, +0, +0, -0.0e+00),
         59.52847891, ""},
        {/*491555*/ DoodsonConstituent(4, +4, -4, +0, +0, +0, +0.0e+00),
         60.00000000, ""},
        {/*493555*/ DoodsonConstituent(4, +4, -2, +0, +0, +0, -0.0e+00),
         60.08213728, ""},
        {/*625655*/ DoodsonConstituent(6, -3, +0, +1, +0, +0, -0.0e+00),
         85.30990488, ""},
        {/*627655*/ DoodsonConstituent(6, -3, +2, +1, +0, +0, -0.0e+00),
         85.39204216, ""},
        {/*635555*/ DoodsonConstituent(6, -2, +0, +0, +0, +0, -0.0e+00),
         85.85427958, ""},
        {/*635755*/ DoodsonConstituent(6, -2, +0, +2, +0, +0, -1.0e+00),
         85.86356325, ""},
        {/*637555*/ DoodsonConstituent(6, -2, +2, +0, +0, +0, -0.0e+00),
         85.93641686, ""},
        {/*645655*/ DoodsonConstituent(6, -1, +0, +1, +0, +0, -0.0e+00),
         86.40793794, ""},
        {/*647455*/ DoodsonConstituent(6, -1, +2, -1, +0, +0, -0.0e+00),
         86.48079155, ""},
        {/*653555*/ DoodsonConstituent(6, +0, -2, +0, +0, +0, -0.0e+00),
         86.87017536, ""},
        {/*655555*/ DoodsonConstituent(6, +0, +0, +0, +0, +0, +0.0e+00),
         86.95231264, ""},
        {/*657555*/ DoodsonConstituent(6, +0, +2, +0, +0, +0, -0.0e+00),
         87.03444992, ""},
        {/*663655*/ DoodsonConstituent(6, +1, -2, +1, +0, +0, -0.0e+00),
         87.42383373, ""},
        {/*665455*/ DoodsonConstituent(6, +1, +0, -1, +0, +0, +1.0e+00),
         87.49668734, ""},
        {/*665655*/ DoodsonConstituent(6, +1, +0, +1, +0, +0, -0.0e+00),
         87.50597101, ""},
        {/*673555*/ DoodsonConstituent(6, +2, -2, +0, +0, +0, -0.0e+00),
         87.96820843, ""},
        {/*675355*/ DoodsonConstituent(6, +2, +0, -2, +0, +0, +1.0e+00),
         88.04106204, ""},
        {/*675555*/ DoodsonConstituent(6, +2, +0, +0, +0, +0, -0.0e+00),
         88.05034571, ""},
        {/*683455*/ DoodsonConstituent(6, +3, -2, -1, +0, +0, -0.0e+00),
         88.51258313, ""},
        {/*685455*/ DoodsonConstituent(6, +3, +0, -1, +0, +0, -0.0e+00),
         88.59472040, ""},
        {/*691555*/ DoodsonConstituent(6, +4, -4, +0, +0, +0, -0.0e+00),
         88.98410421, ""},
        {/*693555*/ DoodsonConstituent(6, +4, -2, +0, +0, +0, -0.0e+00),
         89.06624149, ""},
        {/*835755*/ DoodsonConstituent(8, -2, +0, +2, +0, +0, -0.0e+00),
         114.84766746, ""},
        {/*845655*/ DoodsonConstituent(8, -1, +0, +1, +0, +0, -0.0e+00),
         115.39204216, ""},
        {/*855555*/ DoodsonConstituent(8, +0, +0, +0, +0, +0, +0.0e+00),
         115.93641686, ""},
        {/*863655*/ DoodsonConstituent(8, +1, -2, +1, +0, +0, -0.0e+00),
         116.40793794, ""},
        {/*865455*/ DoodsonConstituent(8, +1, +0, -1, +0, +0, +1.0e+00),
         116.48079155, ""},
        {/*873555*/ DoodsonConstituent(8, +2, -2, +0, +0, +0, -0.0e+00),
         116.95231264, ""},
        {/*875555*/ DoodsonConstituent(8, +2, +0, +0, +0, +0, -0.0e+00),
         117.03444992, ""},
        {/*891555*/ DoodsonConstituent(8, +4, -4, +0, +0, +0, -0.0e+00),
         117.96820843, ""},
        {/*893555*/ DoodsonConstituent(8, +4, -2, +0, +0, +0, -0.0e+00),
         118.05034571, ""},
        {/*895555*/ DoodsonConstituent(8, +4, +0, +0, +0, +0, -0.0e+00),
         118.13248298, ""},
    }}; /* TidalConstituentsArray */
} /* namespace dso */
#endif
