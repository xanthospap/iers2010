/** @file
 * Define a DoodsonConstituent class to represent tidal constituents via
 * Doodson numbers.
 * References:
 * [1]
 * https://ivscc.gsfc.nasa.gov/hfeop_wg/memos/memo-conventions_Ray_2017Dec10.pdf
 */

#ifndef __DOODSON_NUMBER_DEFINES_HPP__
#define __DOODSON_NUMBER_DEFINES_HPP__

#include "datetime/dtcalendar.hpp"
#include "geodesy/units.hpp"
#include <array>
#include <cassert>
#include <cstring>
#include <initializer_list>

namespace dso {

/* @class
 * @ref "Indexing and argument conventions for tides", Richard Ray (GSFC),
 *        2017, available at [1]
 *
 *  Each instance holds the 6 Fundamental Doodson variables, namely:
 *  Table 1: Fundamental variables—Doodson versus Delaunay
 * Variable (k_i)                      Rate (cpd)    Period
 * ---------------------------------------------------------------------------
 * τ  mean lunar time                  9.661 × 10−1  1.03505 d (lunar day)
 * s  mean longitude of moon           3.660 × 10−2  27.32158 d (tropical month)
 * h  mean longitude of sun            2.738 × 10−3  365.2422 d (tropical year)
 * p  mean longitude of lunar perigee  3.095 × 10−4  8.847 y (lunar orbit
 * precession) N' negative longitude of lunar node 1.471 × 10−4  18.61 y
 * (regression of lunar node) ps mean longitude of solar perigee  1.307 × 10−7
 * 21,000 y
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
   *  τ   = ar[0]
   *  s   = ar[1]
   *  h   = ar[2]
   *  p   = ar[3]
   *  N'  = ar[4]
   *  p_s = ar[5]
   */
  explicit DoodsonConstituent(const int *ar = nullptr) {
    if (ar)
      std::memcpy(iar, ar, sizeof(int) * 6);
  }

  /* @brief Constructor from initializer list (of ints)
   *        i.e. enables: DoodsonConstituent d(1,2,3,4,5,6);
   * @param[in] l An initializer list of 6 ints, interpreted as:
   *  τ   = l[0]
   *  s   = l[1]
   *  h   = l[2]
   *  p   = l[3]
   *  N'  = l[4]
   *  p_s = l[5]
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

}; /* DoodsonConstituent */

namespace detail {
struct TidalConstituentsArrayEntry {
  DoodsonConstituent _d;
  const char _n[16];
}; /* TidalConstituentsArrayEntry */
} /* namespace detail */

constexpr static std::array<detail::TidalConstituentsArrayEntry, 64>
    TidalConstituentsArray = {{
        /* long period */
        {DoodsonConstituent(0, 0, 0, 0, 1, 0, 1e0), "om1"},
        {DoodsonConstituent(0, 0, 0, 0, 2, 0, 0e0), "om2"},
        {DoodsonConstituent(0, 0, 1, 0, 0, -1, 0e0), "sa"},
        {DoodsonConstituent(0, 0, 1, 0, 0, 0, 0e0), "sa"},
        {DoodsonConstituent(0, 0, 2, 0, 0, 0, 0e0), "ssa"},
        {DoodsonConstituent(0, 0, 3, 0, 0, -1, 0e0), "sta"},
        {DoodsonConstituent(0, 1, -2, 1, 0, 0, 0e0), "msm"},
        {DoodsonConstituent(0, 1, 0, -1, 0, 0, 0e0), "mm"},
        {DoodsonConstituent(0, 2, -2, 0, 0, 0, 0e0), "msf"},
        {DoodsonConstituent(0, 2, 0, 0, 0, 0, 0e0), "mf"},
        {DoodsonConstituent(0, 3, -2, 1, 0, 0, 0e0), "mstm"},
        {DoodsonConstituent(0, 3, 0, -1, 0, 0, 0e0), "mtm"},
        /* following also goes by the name "msqm" */
        {DoodsonConstituent(0, 4, -2, 0, 0, 0, 0e0), "msq"},
        /* diurnal */
        {DoodsonConstituent(1, -3, 0, 2, 0, 0, .5e0), "2q1"},
        /* following also named "sigma1" */
        {DoodsonConstituent(1, -3, 2, 0, 0, 0, .5e0), "sig1"},
        {DoodsonConstituent(1, -2, 0, 1, 0, 0, .5e0), "q1"},
        /* following also named "rho1" */
        {DoodsonConstituent(1, -2, 2, -1, 0, 0, .5e0), "ro1"},
        {DoodsonConstituent(1, -1, 0, 0, 0, 0, .5e0), "o1"},
        {DoodsonConstituent(1, -1, 2, 0, 0, 0, -.5e0), "tau1"},
        {DoodsonConstituent(1, 0, 0, 1, 0, 0, -.5e0), "m1"},
        {DoodsonConstituent(1, 0, 2, -1, 0, 0, -.5e0), "chi1"},
        {DoodsonConstituent(1, 1, -3, 0, 0, 1, .5e0), "pi1"},
        {DoodsonConstituent(1, 1, -2, 0, 0, 0, .5e0), "p1"},
        {DoodsonConstituent(1, 1, -1, 0, 0, 0, .5e0), "s1"},
        {DoodsonConstituent(1, 1, -1, 0, 0, 1, .5e0), "s1"},
        {DoodsonConstituent(1, 1, 0, 0, 0, 0, -.5e0), "k1"},
        {DoodsonConstituent(1, 1, 1, 0, 0, -1, -.5e0), "psi1"},
        /* following also named "phi1" */
        {DoodsonConstituent(1, 1, 2, 0, 0, 0, -.5e0), "fi1"},
        /* following also named "theta1" */
        {DoodsonConstituent(1, 2, -2, 1, 0, 0, -.5e0), "the1"},
        {DoodsonConstituent(1, 2, 0, -1, 0, 0, -.5e0), "j1"},
        {DoodsonConstituent(1, 3, -2, 0, 0, 0, -.5e0), "so1"},
        {DoodsonConstituent(1, 3, 0, 0, 0, 0, -.5e0), "oo1"},
        {DoodsonConstituent(1, 4, 0, -1, 0, 0, -.5e0), "v1"},
        /* semi-diurnal */
        {DoodsonConstituent(2, -3, 0, 3, 0, 0, 0e0), "3n2"},
        {DoodsonConstituent(2, -3, 2, 1, 0, 0, 0e0), "eps2"},
        {DoodsonConstituent(2, -2, 0, 2, 0, 0, 0e0), "2n2"},
        /* following also named "mi2" */
        {DoodsonConstituent(2, -2, 2, 0, 0, 0, 0e0), "mu2"},
        {DoodsonConstituent(2, -1, 0, 1, 0, 0, 0e0), "n2"},
        /* following also named "ni2" */
        {DoodsonConstituent(2, -1, 2, -1, 0, 0, 0e0), "nu2"},
        {DoodsonConstituent(2, 0, -2, 2, 0, 0, 1e0), "gam2"},
        {DoodsonConstituent(2, 0, -1, 0, 0, 1, 1e0), "alf2"},
        {DoodsonConstituent(2, 0, 0, 0, 0, 0, 0e0), "m2"},
        {DoodsonConstituent(2, 0, 1, 0, 0, -1, 0e0), "bet2"},
        {DoodsonConstituent(2, 0, 2, 0, 0, 0, 0e0), "dlt2"},
        /* following also named "lmb2" and "lambda2" */
        {DoodsonConstituent(2, 1, -2, 1, 0, 0, 1e0), "la2"},
        {DoodsonConstituent(2, 1, 0, -1, 0, 0, 0e0), "l2"},
        {DoodsonConstituent(2, 2, -4, 0, 0, 2, 0e0), "2t2"},
        {DoodsonConstituent(2, 2, -3, 0, 0, 1, 0e0), "t2"},
        {DoodsonConstituent(2, 2, -2, 0, 0, 0, 0e0), "s2"},
        {DoodsonConstituent(2, 2, -1, 0, 0, -1, 1e0), "r2"},
        {DoodsonConstituent(2, 2, 0, 0, 0, 0, 0e0), "k2"},
        {DoodsonConstituent(2, 3, -2, 1, 0, 0, 0e0), "ksi2"},
        {DoodsonConstituent(2, 3, 0, -1, 0, 0, 0e0), "eta2"},
        /* non-linear */
        {DoodsonConstituent(3, 0, 0, 0, 0, 0, 1e0), "m3"},
        {DoodsonConstituent(3, 3, -4, 0, 0, 0, 1e0), "t3"},
        {DoodsonConstituent(3, 3, -3, 0, 0, 0, 1e0), "s3"},
        {DoodsonConstituent(3, 3, -2, 0, 0, 0, 1e0), "r3"},
        {DoodsonConstituent(4, -2, 0, 2, 0, 0, 0e0), "n4"},
        {DoodsonConstituent(4, -1, 0, 1, 0, 0, 0e0), "mn4"},
        {DoodsonConstituent(4, 0, 0, 0, 0, 0, 0e0), "m4"},
        {DoodsonConstituent(4, 2, -2, 0, 0, 0, 0e0), "ms4"},
        {DoodsonConstituent(4, 4, -4, 0, 0, 0, 0e0), "s4"},
        {DoodsonConstituent(6, 0, 0, 0, 0, 0, 0e0), "m6"},
        {DoodsonConstituent(8, 0, 0, 0, 0, 0, 0e0), "m8"}
        /* {"5a0.555", "s5"}, */
        /* {"6bz.555", "s6"}, */
    }}; /* TidalConstituentsArray */

} /* namespace dso */
#endif
