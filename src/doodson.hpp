/** @file
 * Define a DoodsonConstituent class to represent tidal constituents via
 * Doodson numbers.
 * References:
 * [1]
 * https://ivscc.gsfc.nasa.gov/hfeop_wg/memos/memo-conventions_Ray_2017Dec10.pdf
 *
 * [2] Cartwright, D.E. and Edden, A.C. (1973), Corrected Tables of Tidal 
 * Harmonics. Geophysical Journal of the Royal Astronomical Society, 33: 
 * 253-264. https://doi.org/10.1111/j.1365-246X.1973.tb03420.x
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

  double pifactor() const noexcept {return pifac;}

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
  double _hf; /* from IERS 2010, Table 6.7 and [2] */
  const char _n[8]; /* name */
}; /* TidalConstituentsArrayEntry */
} /* namespace detail */

/** Given a tidal constituent name, return its details */
const detail::TidalConstituentsArrayEntry *get_wave(const char *name) noexcept;

/**
 * Refernces: 
 * 1. IERS 2010
 * 2. https://www5.obs-mip.fr/wp-content-omp/uploads/sites/12/2016/10/ttb-1.pdf
 * 3. See also the Section Chapter V. "Harmonic tidal equation" of Bernard 
 * Simon, Coastal Tides, Institut océanographique, Fondation Albert Ier, 
 * Prince de Monaco, 2013
 */
constexpr static std::array<detail::TidalConstituentsArrayEntry, 384>
    TidalConstituentsArray = {{
{/*055555*/DoodsonConstituent( 0,+0,+0,+0,+0,+0, +0.0e+00), -999.00000000, -0.31455, ""},
{/*055565*/DoodsonConstituent( 0,+0,+0,+0,+1,+0, +1.0e+00),   0.00221000, +0.02793, "Om1"},
{/*055575*/DoodsonConstituent( 0,+0,+0,+0,+2,+0, +0.0e+00),   0.00441000, -0.00028, "Om2"},
{/*055765*/DoodsonConstituent( 0,+0,+0,+2,+1,+0, +1.0e+00), -999.00000000, +0.00004, ""},
{/*056544*/DoodsonConstituent( 0,+0,+1,+0,-1,-1, +0.0e+00), -999.00000000, -0.00004, ""},
{/*056554*/DoodsonConstituent( 0,+0,+1,+0,+0,-1, +0.0e+00),   0.04107000, -0.00492, "Sa"},
//{/*056555*/DoodsonConstituent( 0,+0,+1,+0,+0,+0, +0.0e+00),   0.04106864, -999.00000, "sa"},
{/*056556*/DoodsonConstituent( 0,+0,+1,+0,+0,+1, +1.0e+00), -999.00000000, +0.00026, ""},
{/*056564*/DoodsonConstituent( 0,+0,+1,+0,+1,-1, +1.0e+00), -999.00000000, +0.00005, ""},
{/*057345*/DoodsonConstituent( 0,+0,+2,-2,-1,+0, +1.0e+00), -999.00000000, +0.00002, ""},
{/*057355*/DoodsonConstituent( 0,+0,+2,-2,+0,+0, +0.0e+00), -999.00000000, -0.00032, ""},
{/*057553*/DoodsonConstituent( 0,+0,+2,+0,+0,-2, +0.0e+00), -999.00000000, -0.00012, ""},
{/*057555*/DoodsonConstituent( 0,+0,+2,+0,+0,+0, +0.0e+00),   0.08213728, -0.03100, "Ssa"},
{/*057565*/DoodsonConstituent( 0,+0,+2,+0,+1,+0, +1.0e+00),   0.08434000, +0.00077, ""},
{/*057575*/DoodsonConstituent( 0,+0,+2,+0,+2,+0, +1.0e+00), -999.00000000, +0.00017, ""},
{/*058554*/DoodsonConstituent( 0,+0,+3,+0,+0,-1, +0.0e+00),   0.12320000, -0.00181, "Sta"},
{/*058564*/DoodsonConstituent( 0,+0,+3,+0,+1,-1, +1.0e+00), -999.00000000, +0.00003, ""},
{/*059553*/DoodsonConstituent( 0,+0,+4,+0,+0,-2, +0.0e+00), -999.00000000, -0.00007, ""},
{/*062646*/DoodsonConstituent( 0,+1,-3,+1,-1,+1, +1.0e+00), -999.00000000, +0.00002, ""},
{/*062656*/DoodsonConstituent( 0,+1,-3,+1,+0,+1, +0.0e+00), -999.00000000, -0.00029, ""},
{/*063435*/DoodsonConstituent( 0,+1,-2,-1,-2,+0, +1.0e+00), -999.00000000, +0.00002, ""},
{/*063445*/DoodsonConstituent( 0,+1,-2,-1,-1,+0, +1.0e+00), -999.00000000, +0.00007, ""},
{/*063645*/DoodsonConstituent( 0,+1,-2,+1,-1,+0, +1.0e+00), -999.00000000, +0.00048, ""},
{/*063655*/DoodsonConstituent( 0,+1,-2,+1,+0,+0, +0.0e+00),   0.47152109, -0.00673, "Msm"},
{/*063665*/DoodsonConstituent( 0,+1,-2,+1,+1,+0, +1.0e+00), -999.00000000, +0.00044, ""},
{/*064456*/DoodsonConstituent( 0,+1,-1,-1,+0,+1, +0.0e+00), -999.00000000, -0.00022, ""},
{/*064555*/DoodsonConstituent( 0,+1,-1,+0,+0,+0, +1.0e+00), -999.00000000, +0.00020, ""},
{/*064654*/DoodsonConstituent( 0,+1,-1,+1,+0,-1, +1.0e+00), -999.00000000, +0.00005, ""},
{/*065445*/DoodsonConstituent( 0,+1,+0,-1,-1,+0, +1.0e+00),   0.54217000, +0.00231, ""},
{/*065455*/DoodsonConstituent( 0,+1,+0,-1,+0,+0, +0.0e+00),   0.54437470, -0.03518, "Mm"},
{/*065465*/DoodsonConstituent( 0,+1,+0,-1,+1,+0, +1.0e+00),   0.54658000, +0.00229, ""},
{/*065555*/DoodsonConstituent( 0,+1,+0,+0,+0,+0, +0.0e+00), -999.00000000, -0.00375, ""},
{/*065655*/DoodsonConstituent( 0,+1,+0,+1,+0,+0, +1.0e+00),   0.55366000, +0.00188, ""},
{/*065665*/DoodsonConstituent( 0,+1,+0,+1,+1,+0, +1.0e+00), -999.00000000, +0.00077, ""},
{/*065675*/DoodsonConstituent( 0,+1,+0,+1,+2,+0, +1.0e+00), -999.00000000, +0.00021, ""},
{/*066454*/DoodsonConstituent( 0,+1,+1,-1,+0,-1, +1.0e+00), -999.00000000, +0.00018, ""},
{/*067455*/DoodsonConstituent( 0,+1,+2,-1,+0,+0, +1.0e+00), -999.00000000, +0.00049, ""},
{/*067465*/DoodsonConstituent( 0,+1,+2,-1,+1,+0, +1.0e+00), -999.00000000, +0.00024, ""},
{/*067475*/DoodsonConstituent( 0,+1,+2,-1,+2,+0, +1.0e+00), -999.00000000, +0.00004, ""},
{/*068454*/DoodsonConstituent( 0,+1,+3,-1,+0,-1, +1.0e+00), -999.00000000, +0.00002, ""},
{/*071755*/DoodsonConstituent( 0,+2,-4,+2,+0,+0, +0.0e+00), -999.00000000, -0.00011, ""},
{/*072556*/DoodsonConstituent( 0,+2,-3,+0,+0,+1, +0.0e+00), -999.00000000, -0.00038, ""},
{/*072566*/DoodsonConstituent( 0,+2,-3,+0,+1,+1, +1.0e+00), -999.00000000, +0.00002, ""},
{/*073545*/DoodsonConstituent( 0,+2,-2,+0,-1,+0, +0.0e+00), -999.00000000, -0.00042, ""},
{/*073555*/DoodsonConstituent( 0,+2,-2,+0,+0,+0, +0.0e+00),   1.01589579, -0.00583, "Msf"},
{/*073565*/DoodsonConstituent( 0,+2,-2,+0,+1,+0, +1.0e+00), -999.00000000, +0.00038, ""},
{/*073755*/DoodsonConstituent( 0,+2,-2,+2,+0,+0, +1.0e+00), -999.00000000, +0.00004, ""},
{/*074356*/DoodsonConstituent( 0,+2,-1,-2,+0,+1, +0.0e+00), -999.00000000, -0.00004, ""},
{/*074455*/DoodsonConstituent( 0,+2,-1,-1,+0,+0, +1.0e+00), -999.00000000, +0.00003, ""},
{/*074554*/DoodsonConstituent( 0,+2,-1,+0,+0,-1, +1.0e+00), -999.00000000, +0.00006, ""},
{/*074556*/DoodsonConstituent( 0,+2,-1,+0,+0,+1, +0.0e+00), -999.00000000, -0.00020, ""},
{/*074566*/DoodsonConstituent( 0,+2,-1,+0,+1,+1, +0.0e+00), -999.00000000, -0.00004, ""},
{/*075345*/DoodsonConstituent( 0,+2,+0,-2,-1,+0, +1.0e+00), -999.00000000, +0.00015, ""},
{/*075355*/DoodsonConstituent( 0,+2,+0,-2,+0,+0, +0.0e+00),   1.08875000, -0.00288, ""},
{/*075365*/DoodsonConstituent( 0,+2,+0,-2,+1,+0, +1.0e+00), -999.00000000, +0.00019, ""},
{/*075555*/DoodsonConstituent( 0,+2,+0,+0,+0,+0, +0.0e+00),   1.09803306, -0.06663, "Mf"},
{/*075565*/DoodsonConstituent( 0,+2,+0,+0,+1,+0, +0.0e+00),   1.10024000, -0.02762, ""},
{/*075575*/DoodsonConstituent( 0,+2,+0,+0,+2,+0, +0.0e+00),   1.10245000, -0.00258, ""},
{/*075585*/DoodsonConstituent( 0,+2,+0,+0,+3,+0, +1.0e+00), -999.00000000, +0.00006, ""},
{/*076354*/DoodsonConstituent( 0,+2,+1,-2,+0,-1, +1.0e+00), -999.00000000, +0.00003, ""},
{/*076554*/DoodsonConstituent( 0,+2,+1,+0,+0,-1, +1.0e+00), -999.00000000, +0.00023, ""},
{/*076564*/DoodsonConstituent( 0,+2,+1,+0,+1,-1, +1.0e+00), -999.00000000, +0.00006, ""},
{/*077355*/DoodsonConstituent( 0,+2,+2,-2,+0,+0, +1.0e+00), -999.00000000, +0.00020, ""},
{/*077365*/DoodsonConstituent( 0,+2,+2,-2,+1,+0, +1.0e+00), -999.00000000, +0.00008, ""},
{/*077575*/DoodsonConstituent( 0,+2,+2,+0,+2,+0, +1.0e+00), -999.00000000, +0.00003, ""},
{/*081655*/DoodsonConstituent( 0,+3,-4,+1,+0,+0, +0.0e+00), -999.00000000, -0.00017, ""},
{/*082456*/DoodsonConstituent( 0,+3,-3,-1,+0,+1, +0.0e+00), -999.00000000, -0.00007, ""},
{/*082656*/DoodsonConstituent( 0,+3,-3,+1,+0,+1, +0.0e+00), -999.00000000, -0.00011, ""},
{/*082666*/DoodsonConstituent( 0,+3,-3,+1,+1,+1, +0.0e+00), -999.00000000, -0.00004, ""},
{/*083445*/DoodsonConstituent( 0,+3,-2,-1,-1,+0, +0.0e+00), -999.00000000, -0.00009, ""},
{/*083455*/DoodsonConstituent( 0,+3,-2,-1,+0,+0, +0.0e+00), -999.00000000, -0.00092, ""},
{/*083465*/DoodsonConstituent( 0,+3,-2,-1,+1,+0, +1.0e+00), -999.00000000, +0.00006, ""},
{/*083655*/DoodsonConstituent( 0,+3,-2,+1,+0,+0, +0.0e+00),   1.56955415, -0.00242, "Mstm"},
{/*083665*/DoodsonConstituent( 0,+3,-2,+1,+1,+0, +0.0e+00), -999.00000000, -0.00100, ""},
{/*083675*/DoodsonConstituent( 0,+3,-2,+1,+2,+0, +0.0e+00), -999.00000000, -0.00009, ""},
{/*084456*/DoodsonConstituent( 0,+3,-1,-1,+0,+1, +0.0e+00), -999.00000000, -0.00013, ""},
{/*084466*/DoodsonConstituent( 0,+3,-1,-1,+1,+1, +0.0e+00), -999.00000000, -0.00004, ""},
{/*084565*/DoodsonConstituent( 0,+3,-1,+0,+1,+0, +1.0e+00), -999.00000000, +0.00003, ""},
{/*084654*/DoodsonConstituent( 0,+3,-1,+1,+0,-1, +1.0e+00), -999.00000000, +0.00003, ""},
{/*085255*/DoodsonConstituent( 0,+3,+0,-3,+0,+0, +0.0e+00), -999.00000000, -0.00023, ""},
{/*085264*/DoodsonConstituent( 0,+3,+0,-3,+1,-1, +1.0e+00), -999.00000000, +0.00004, ""},
{/*085266*/DoodsonConstituent( 0,+3,+0,-3,+1,+1, +1.0e+00), -999.00000000, +0.00004, ""},
{/*085455*/DoodsonConstituent( 0,+3,+0,-1,+0,+0, +0.0e+00),   1.64240776, -0.01276, "Mtm"},
{/*085465*/DoodsonConstituent( 0,+3,+0,-1,+1,+0, +0.0e+00),   1.64462000, -0.00529, ""},
{/*085475*/DoodsonConstituent( 0,+3,+0,-1,+2,+0, +0.0e+00), -999.00000000, -0.00051, ""},
{/*085675*/DoodsonConstituent( 0,+3,+0,+1,+2,+0, +1.0e+00), -999.00000000, +0.00005, ""},
{/*085685*/DoodsonConstituent( 0,+3,+0,+1,+3,+0, +1.0e+00), -999.00000000, +0.00002, ""},
{/*086454*/DoodsonConstituent( 0,+3,+1,-1,+0,-1, +1.0e+00), -999.00000000, +0.00011, ""},
{/*091555*/DoodsonConstituent( 0,+4,-4,+0,+0,+0, +0.0e+00), -999.00000000, -0.00009, ""},
{/*091755*/DoodsonConstituent( 0,+4,-4,+2,+0,+0, +0.0e+00), -999.00000000, -0.00006, ""},
{/*091765*/DoodsonConstituent( 0,+4,-4,+2,+1,+0, +0.0e+00), -999.00000000, -0.00003, ""},
{/*092556*/DoodsonConstituent( 0,+4,-3,+0,+0,+1, +0.0e+00), -999.00000000, -0.00014, ""},
{/*092566*/DoodsonConstituent( 0,+4,-3,+0,+1,+1, +0.0e+00), -999.00000000, -0.00006, ""},
{/*093355*/DoodsonConstituent( 0,+4,-2,-2,+0,+0, +0.0e+00), -999.00000000, -0.00011, ""},
{/*093555*/DoodsonConstituent( 0,+4,-2,+0,+0,+0, +0.0e+00),   2.11394000, -0.00204, "Msqm"},
{/*093565*/DoodsonConstituent( 0,+4,-2,+0,+1,+0, +0.0e+00), -999.00000000, -0.00084, ""},
{/*093575*/DoodsonConstituent( 0,+4,-2,+0,+2,+0, +0.0e+00), -999.00000000, -0.00008, ""},
{/*093655*/DoodsonConstituent( 0,+4,-2,+1,+0,+0, +0.0e+00), -999.00000000, -0.00242, ""},
{/*094356*/DoodsonConstituent( 0,+4,-1,-2,+0,+1, +0.0e+00), -999.00000000, -0.00003, ""},
{/*094554*/DoodsonConstituent( 0,+4,-1,+0,+0,-1, +1.0e+00), -999.00000000, +0.00003, ""},
{/*094555*/DoodsonConstituent( 0,+4,-1,+0,+0,+0, +1.0e+00), -999.00000000, +0.00007, ""},
{/*095355*/DoodsonConstituent( 0,+4,+0,-2,+0,+0, +0.0e+00),   2.18678246, -0.00169, ""},
{/*095365*/DoodsonConstituent( 0,+4,+0,-2,+1,+0, +0.0e+00), -999.00000000, -0.00070, ""},
{/*095375*/DoodsonConstituent( 0,+4,+0,-2,+2,+0, +0.0e+00), -999.00000000, -0.00007, ""},
//{/*105555*/DoodsonConstituent( 1,-5,+0,+0,+0,+0, +2.5e+00),  11.74696945, -999.00000, ""},
{/*115845*/DoodsonConstituent( 1,-4,+0,+3,-1,+0, -5.0e-01), -999.00000000, -0.00014, ""},
{/*115855*/DoodsonConstituent( 1,-4,+0,+3,+0,+0, -5.0e-01), -999.00000000, -0.00075, ""},
{/*116656*/DoodsonConstituent( 1,-4,+1,+1,+0,+1, +5.0e-01), -999.00000000, +0.00004, ""},
{/*117645*/DoodsonConstituent( 1,-4,+2,+1,-1,+0, -5.0e-01), -999.00000000, -0.00037, ""},
{/*117655*/DoodsonConstituent( 1,-4,+2,+1,+0,+0, -5.0e-01), -999.00000000, -0.00194, ""},
{/*118654*/DoodsonConstituent( 1,-4,+3,+1,+0,-1, -5.0e-01), -999.00000000, -0.00015, ""},
{/*119445*/DoodsonConstituent( 1,-4,+4,-1,-1,+0, -5.0e-01), -999.00000000, -0.00007, ""},
{/*119455*/DoodsonConstituent( 1,-4,+4,-1,+0,+0, -5.0e-01), -999.00000000, -0.00037, ""},
{/*124756*/DoodsonConstituent( 1,-3,-1,+2,+0,+1, +5.0e-01), -999.00000000, +0.00009, ""},
{/*125535*/DoodsonConstituent( 1,-3,+0,+0,-2,+0, +5.0e-01), -999.00000000, +0.00004, ""},
{/*125735*/DoodsonConstituent( 1,-3,+0,+2,-2,+0, +5.0e-01), -999.00000000, +0.00003, ""},
{/*125745*/DoodsonConstituent( 1,-3,+0,+2,-1,+0, -5.0e-01), -999.00000000, -0.00125, ""},
{/*125755*/DoodsonConstituent( 1,-3,+0,+2,+0,+0, -5.0e-01),  12.85428618, -0.00664, "2Q1"},
{/*126556*/DoodsonConstituent( 1,-3,+1,+0,+0,+1, +5.0e-01), -999.00000000, +0.00011, ""},
{/*126655*/DoodsonConstituent( 1,-3,+1,+1,+0,+0, +5.0e-01), -999.00000000, +0.00007, ""},
{/*126754*/DoodsonConstituent( 1,-3,+1,+2,+0,-1, -5.0e-01), -999.00000000, -0.00010, ""},
{/*127535*/DoodsonConstituent( 1,-3,+2,+0,-2,+0, +5.0e-01), -999.00000000, +0.00004, ""},
{/*127545*/DoodsonConstituent( 1,-3,+2,+0,-1,+0, -5.0e-01), -999.00000000, -0.00151, ""},
{/*127555*/DoodsonConstituent( 1,-3,+2,+0,+0,+0, -5.0e-01),  12.92713979, -0.00802, "sig1"},
{/*127755*/DoodsonConstituent( 1,-3,+2,+2,+0,+0, +5.0e-01), -999.00000000, +0.00007, ""},
{/*128544*/DoodsonConstituent( 1,-3,+3,+0,-1,-1, -5.0e-01), -999.00000000, -0.00010, ""},
{/*128554*/DoodsonConstituent( 1,-3,+3,+0,+0,-1, -5.0e-01), -999.00000000, -0.00054, ""},
{/*129345*/DoodsonConstituent( 1,-3,+4,-2,-1,+0, -5.0e-01), -999.00000000, -0.00005, ""},
{/*129355*/DoodsonConstituent( 1,-3,+4,-2,+0,+0, -5.0e-01), -999.00000000, -0.00024, ""},
{/*129555*/DoodsonConstituent( 1,-3,+4,+0,+0,+0, +5.0e-01), -999.00000000, +0.00008, ""},
{/*129565*/DoodsonConstituent( 1,-3,+4,+0,+1,+0, -5.0e-01), -999.00000000, -0.00003, ""},
{/*133635*/DoodsonConstituent( 1,-2,-2,+1,-2,+0, +5.0e-01), -999.00000000, +0.00004, ""},
{/*133855*/DoodsonConstituent( 1,-2,-2,+3,+0,+0, +5.0e-01), -999.00000000, +0.00016, ""},
{/*134646*/DoodsonConstituent( 1,-2,-1,+1,-1,+1, +5.0e-01), -999.00000000, +0.00007, ""},
{/*135425*/DoodsonConstituent( 1,-2,+0,-1,-3,+0, +5.0e-01), -999.00000000, +0.00004, ""},
{/*135435*/DoodsonConstituent( 1,-2,+0,-1,-2,+0, +5.0e-01), -999.00000000, +0.00019, ""},
{/*135556*/DoodsonConstituent( 1,-2,+0,+0,+0,+1, -5.0e-01), -999.00000000, -0.00004, ""},
{/*135635*/DoodsonConstituent( 1,-2,+0,+1,-2,+0, +5.0e-01), -999.00000000, +0.00029, ""},
{/*135645*/DoodsonConstituent( 1,-2,+0,+1,-1,+0, -5.0e-01),  13.39645000, -0.00947, "sig1"},
{/*135655*/DoodsonConstituent( 1,-2,+0,+1,+0,+0, -5.0e-01),  13.39866088, -0.05020, "Q1"},
{/*135855*/DoodsonConstituent( 1,-2,+0,+3,+0,+0, +5.0e-01), -999.00000000, +0.00014, ""},
{/*136456*/DoodsonConstituent( 1,-2,+1,-1,+0,+1, +5.0e-01), -999.00000000, +0.00009, ""},
{/*136545*/DoodsonConstituent( 1,-2,+1,+0,-1,+0, +5.0e-01), -999.00000000, +0.00005, ""},
{/*136555*/DoodsonConstituent( 1,-2,+1,+0,+0,+0, +5.0e-01), -999.00000000, +0.00027, ""},
{/*136644*/DoodsonConstituent( 1,-2,+1,+1,-1,-1, -5.0e-01), -999.00000000, -0.00008, ""},
{/*136654*/DoodsonConstituent( 1,-2,+1,+1,+0,-1, -5.0e-01), -999.00000000, -0.00046, ""},
{/*137435*/DoodsonConstituent( 1,-2,+2,-1,-2,+0, +5.0e-01), -999.00000000, +0.00005, ""},
{/*137442*/DoodsonConstituent( 1,-2,+2,-1,-1,-3, -5.0e-01), -999.00000000, -0.00180, ""},
{/*137445*/DoodsonConstituent( 1,-2,+2,-1,-1,+0, -5.0e-01), -999.00000000, -0.00180, ""},
{/*137455*/DoodsonConstituent( 1,-2,+2,-1,+0,+0, -5.0e-01),  13.47151449, -0.00954, "rho1"},
{/*137655*/DoodsonConstituent( 1,-2,+2,+1,+0,+0, +5.0e-01), -999.00000000, +0.00055, ""},
{/*137665*/DoodsonConstituent( 1,-2,+2,+1,+1,+0, -5.0e-01), -999.00000000, -0.00017, ""},
{/*138444*/DoodsonConstituent( 1,-2,+3,-1,-1,-1, -5.0e-01), -999.00000000, -0.00008, ""},
{/*138454*/DoodsonConstituent( 1,-2,+3,-1,+0,-1, -5.0e-01), -999.00000000, -0.00044, ""},
{/*138654*/DoodsonConstituent( 1,-2,+3,+1,+0,-1, +5.0e-01), -999.00000000, +0.00004, ""},
{/*139455*/DoodsonConstituent( 1,-2,+4,-1,+0,+0, +5.0e-01), -999.00000000, +0.00012, ""},
{/*143535*/DoodsonConstituent( 1,-1,-2,+0,-2,+0, +5.0e-01), -999.00000000, +0.00011, ""},
{/*143745*/DoodsonConstituent( 1,-1,-2,+2,-1,+0, +5.0e-01), -999.00000000, +0.00014, ""},
{/*143755*/DoodsonConstituent( 1,-1,-2,+2,+0,+0, +5.0e-01), -999.00000000, +0.00079, ""},
{/*144546*/DoodsonConstituent( 1,-1,-1,+0,-1,+1, +5.0e-01), -999.00000000, +0.00011, ""},
{/*144556*/DoodsonConstituent( 1,-1,-1,+0,+0,+1, +5.0e-01), -999.00000000, +0.00090, ""},
{/*144655*/DoodsonConstituent( 1,-1,-1,+1,+0,+0, -5.0e-01), -999.00000000, -0.00004, ""},
{/*145535*/DoodsonConstituent( 1,-1,+0,+0,-2,+0, +5.0e-01), -999.00000000, +0.00152, ""},
{/*145545*/DoodsonConstituent( 1,-1,+0,+0,-1,+0, -5.0e-01),  13.94083000, -0.04945, ""},
{/*145555*/DoodsonConstituent( 1,-1,+0,+0,+0,+0, -5.0e-01),  13.94303557, -0.26221, "O1"},
{/*145745*/DoodsonConstituent( 1,-1,+0,+2,-1,+0, -5.0e-01), -999.00000000, -0.00005, ""},
{/*145755*/DoodsonConstituent( 1,-1,+0,+2,+0,+0, +5.0e-01), -999.00000000, +0.00170, ""},
{/*145765*/DoodsonConstituent( 1,-1,+0,+2,+1,+0, +5.0e-01), -999.00000000, +0.00028, ""},
{/*146544*/DoodsonConstituent( 1,-1,+1,+0,-1,-1, -5.0e-01), -999.00000000, -0.00008, ""},
{/*146554*/DoodsonConstituent( 1,-1,+1,+0,+0,-1, -5.0e-01), -999.00000000, -0.00076, ""},
//{/*146555*/DoodsonConstituent( 1,-1,+1,+0,+0,+0, -0.0e+00),  13.98410421, -999.00000, ""},
{/*147355*/DoodsonConstituent( 1,-1,+2,-2,+0,+0, +5.0e-01), -999.00000000, +0.00015, ""},
{/*147545*/DoodsonConstituent( 1,-1,+2,+0,-1,+0, -5.0e-01), -999.00000000, -0.00010, ""},
{/*147555*/DoodsonConstituent( 1,-1,+2,+0,+0,+0, +5.0e-01),  14.02517285, +0.00343, "tau1"},
{/*147565*/DoodsonConstituent( 1,-1,+2,+0,+1,+0, -5.0e-01), -999.00000000, -0.00075, ""},
{/*147575*/DoodsonConstituent( 1,-1,+2,+0,+2,+0, -5.0e-01), -999.00000000, -0.00005, ""},
{/*148554*/DoodsonConstituent( 1,-1,+3,+0,+0,-1, +5.0e-01), -999.00000000, +0.00023, ""},
{/*149355*/DoodsonConstituent( 1,-1,+4,-2,+0,+0, +5.0e-01), -999.00000000, +0.00006, ""},
{/*152656*/DoodsonConstituent( 1,+0,-3,+1,+0,+1, +5.0e-01), -999.00000000, +0.00009, ""},
{/*153645*/DoodsonConstituent( 1,+0,-2,+1,-1,+0, +5.0e-01), -999.00000000, +0.00044, ""},
{/*153655*/DoodsonConstituent( 1,+0,-2,+1,+0,+0, +5.0e-01),  14.41456000, +0.00194, ""},
{/*154555*/DoodsonConstituent( 1,+0,-1,+0,+0,+0, -5.0e-01), -999.00000000, -0.00004, ""},
{/*154656*/DoodsonConstituent( 1,+0,-1,+1,+0,+1, -5.0e-01), -999.00000000, -0.00010, ""},
{/*155435*/DoodsonConstituent( 1,+0,+0,-1,-2,+0, -5.0e-01), -999.00000000, -0.00012, ""},
{/*155445*/DoodsonConstituent( 1,+0,+0,-1,-1,+0, +5.0e-01),  14.48520000, +0.00137, ""},
{/*155455*/DoodsonConstituent( 1,+0,+0,-1,+0,+0, +5.0e-01),  14.48741000, +0.00741, ""},
{/*155555*/DoodsonConstituent( 1,+0,+0,+0,+0,+0, -5.0e-01), -999.00000000, -0.00399, ""},
{/*155645*/DoodsonConstituent( 1,+0,+0,+1,-1,+0, -5.0e-01), -999.00000000, -0.00059, ""},
{/*155655*/DoodsonConstituent( 1,+0,+0,+1,+0,+0, +5.0e-01),  14.49669394, +0.02062, "M1"},
{/*155665*/DoodsonConstituent( 1,+0,+0,+1,+1,+0, +5.0e-01),  14.49890000, +0.00414, ""},
{/*155675*/DoodsonConstituent( 1,+0,+0,+1,+2,+0, -5.0e-01), -999.00000000, -0.00011, ""},
{/*156555*/DoodsonConstituent( 1,+0,+1,+0,+0,+0, -5.0e-01), -999.00000000, -0.00012, ""},
{/*156654*/DoodsonConstituent( 1,+0,+1,+1,+0,-1, +5.0e-01), -999.00000000, +0.00013, ""},
{/*157445*/DoodsonConstituent( 1,+0,+2,-1,-1,+0, -5.0e-01), -999.00000000, -0.00011, ""},
{/*157455*/DoodsonConstituent( 1,+0,+2,-1,+0,+0, +5.0e-01),  14.56954755, +0.00394, "chi1"},
{/*157465*/DoodsonConstituent( 1,+0,+2,-1,+1,+0, +5.0e-01),  14.57176000, +0.00087, ""},
{/*158464*/DoodsonConstituent( 1,+0,+3,-1,+1,-1, +5.0e-01), -999.00000000, +0.00004, ""},
{/*159454*/DoodsonConstituent( 1,+0,+4,-1,+0,-1, +5.0e-01), -999.00000000, +0.00017, ""},
{/*161557*/DoodsonConstituent( 1,+1,-4,+0,+0,+2, -5.0e-01), -999.00000000, -0.00029, ""},
{/*162546*/DoodsonConstituent( 1,+1,-3,+0,-1,+1, +5.0e-01), -999.00000000, +0.00006, ""},
{/*162556*/DoodsonConstituent( 1,+1,-3,+0,+0,+1, -5.0e-01),  14.91786468, -0.00714, "pi1"},
{/*163535*/DoodsonConstituent( 1,+1,-2,+0,-2,+0, -5.0e-01), -999.00000000, -0.00010, ""},
{/*163545*/DoodsonConstituent( 1,+1,-2,+0,-1,+0, +5.0e-01),  14.95673000, +0.00137, ""},
{/*163555*/DoodsonConstituent( 1,+1,-2,+0,+0,+0, -5.0e-01),  14.95893136, -0.12203, "P1"},
{/*163557*/DoodsonConstituent( 1,+1,-2,+0,+0,+2, +5.0e-01), -999.00000000, +0.00005, ""},
{/*163755*/DoodsonConstituent( 1,+1,-2,+2,+0,+0, +5.0e-01), -999.00000000, +0.00018, ""},
{/*163765*/DoodsonConstituent( 1,+1,-2,+2,+1,+0, +5.0e-01), -999.00000000, +0.00004, ""},
{/*164554*/DoodsonConstituent( 1,+1,-1,+0,+0,-1, +5.0e-01),  15.00000000, +0.00102, ""},
//{/*164555*/DoodsonConstituent( 1,+1,-1,+0,+0,+0, +1.0e+00),  15.00000000, -999.00000, "s1"},
{/*164556*/DoodsonConstituent( 1,+1,-1,+0,+0,+1, +5.0e-01),  15.00000000, +0.00289, "S1"},
{/*164566*/DoodsonConstituent( 1,+1,-1,+0,+1,+1, -5.0e-01), -999.00000000, -0.00008, ""},
{/*165345*/DoodsonConstituent( 1,+1,+0,-2,-1,+0, +5.0e-01),  15.02958000, +0.00007, ""},
{/*165535*/DoodsonConstituent( 1,+1,+0,+0,-2,+0, +5.0e-01),  15.03665000, +0.00005, ""},
{/*165545*/DoodsonConstituent( 1,+1,+0,+0,-1,+0, -5.0e-01),  15.03886000, -0.00730, "K1-"},
{/*165555*/DoodsonConstituent( 1,+1,+0,+0,+0,+0, +5.0e-01),  15.04106864, +0.36878, "K1"},
{/*165565*/DoodsonConstituent( 1,+1,+0,+0,+1,+0, +5.0e-01),  15.04328000, +0.05001, "K1+"},
{/*165575*/DoodsonConstituent( 1,+1,+0,+0,+2,+0, -5.0e-01),  15.04548000, -0.00108, ""},
//{/*166455*/DoodsonConstituent( 1,+1,+1,-1,+0,+0, -1.0e+03),  15.07749000, -999.00000, ""},
//{/*166544*/DoodsonConstituent( 1,+1,+1,+0,-1,-1, -1.0e+03),  15.07993000, -999.00000, ""},
{/*166554*/DoodsonConstituent( 1,+1,+1,+0,+0,-1, +5.0e-01),  15.08213532, +0.00293, "psi1"},
//{/*166556*/DoodsonConstituent( 1,+1,+1,+0,+0,+1, -1.0e+03),  15.08214000, -999.00000, ""},
{/*166564*/DoodsonConstituent( 1,+1,+1,+0,+1,-1, +5.0e-01),  15.08434000, +0.00005, ""},
{/*167355*/DoodsonConstituent( 1,+1,+2,-2,+0,+0, +5.0e-01),  15.11392000, +0.00018, ""},
{/*167365*/DoodsonConstituent( 1,+1,+2,-2,+1,+0, +5.0e-01),  15.11613000, +0.00005, ""},
{/*167553*/DoodsonConstituent( 1,+1,+2,+0,+0,-2, +5.0e-01), -999.00000000, +0.00007, ""},
{/*167555*/DoodsonConstituent( 1,+1,+2,+0,+0,+0, +5.0e-01),  15.12320592, +0.00525, "phi1"},
{/*167565*/DoodsonConstituent( 1,+1,+2,+0,+1,+0, -5.0e-01),  15.12542000, -0.00020, ""},
{/*167575*/DoodsonConstituent( 1,+1,+2,+0,+2,+0, -5.0e-01), -999.00000000, -0.00010, ""},
{/*168554*/DoodsonConstituent( 1,+1,+3,+0,+0,-1, +5.0e-01),  15.16427000, +0.00031, ""},
{/*172656*/DoodsonConstituent( 1,+2,-3,+1,+0,+1, +5.0e-01), -999.00000000, +0.00017, ""},
{/*173445*/DoodsonConstituent( 1,+2,-2,-1,-1,+0, +5.0e-01), -999.00000000, +0.00012, ""},
{/*173645*/DoodsonConstituent( 1,+2,-2,+1,-1,+0, -5.0e-01), -999.00000000, -0.00012, ""},
{/*173655*/DoodsonConstituent( 1,+2,-2,+1,+0,+0, +5.0e-01),  15.51258973, +0.00395, "theta1"},
{/*173665*/DoodsonConstituent( 1,+2,-2,+1,+1,+0, +5.0e-01),  15.51480000, +0.00073, ""},
{/*174456*/DoodsonConstituent( 1,+2,-1,-1,+0,+1, +5.0e-01), -999.00000000, +0.00012, ""},
{/*174555*/DoodsonConstituent( 1,+2,-1,+0,+0,+0, -5.0e-01), -999.00000000, -0.00012, ""},
{/*175445*/DoodsonConstituent( 1,+2,+0,-1,-1,+0, -5.0e-01),  15.58323000, -0.00060, ""},
{/*175455*/DoodsonConstituent( 1,+2,+0,-1,+0,+0, +5.0e-01),  15.58544334, +0.02062, "J1"},
{/*175465*/DoodsonConstituent( 1,+2,+0,-1,+1,+0, +5.0e-01),  15.58765000, +0.00409, ""},
{/*175475*/DoodsonConstituent( 1,+2,+0,-1,+2,+0, -5.0e-01), -999.00000000, -0.00007, ""},
{/*175655*/DoodsonConstituent( 1,+2,+0,+1,+0,+0, -5.0e-01), -999.00000000, -0.00032, ""},
{/*175665*/DoodsonConstituent( 1,+2,+0,+1,+1,+0, -5.0e-01), -999.00000000, -0.00020, ""},
{/*175675*/DoodsonConstituent( 1,+2,+0,+1,+2,+0, -5.0e-01), -999.00000000, -0.00012, ""},
{/*176454*/DoodsonConstituent( 1,+2,+1,-1,+0,-1, -5.0e-01), -999.00000000, -0.00011, ""},
{/*177455*/DoodsonConstituent( 1,+2,+2,-1,+0,+0, -5.0e-01), -999.00000000, -0.00008, ""},
{/*177465*/DoodsonConstituent( 1,+2,+2,-1,+1,+0, -5.0e-01), -999.00000000, -0.00006, ""},
{/*181755*/DoodsonConstituent( 1,+3,-4,+2,+0,+0, +5.0e-01), -999.00000000, +0.00006, ""},
{/*182556*/DoodsonConstituent( 1,+3,-3,+0,+0,+1, +5.0e-01), -999.00000000, +0.00023, ""},
{/*182566*/DoodsonConstituent( 1,+3,-3,+0,+1,+1, +5.0e-01), -999.00000000, +0.00004, ""},
{/*183545*/DoodsonConstituent( 1,+3,-2,+0,-1,+0, +5.0e-01), -999.00000000, +0.00011, ""},
{/*183555*/DoodsonConstituent( 1,+3,-2,+0,+0,+0, +5.0e-01),  16.05696443, +0.00342, "So1"},
{/*183565*/DoodsonConstituent( 1,+3,-2,+0,+1,+0, +5.0e-01), -999.00000000, +0.00067, ""},
{/*184554*/DoodsonConstituent( 1,+3,-1,+0,+0,-1, -5.0e-01), -999.00000000, -0.00007, ""},
{/*185345*/DoodsonConstituent( 1,+3,+0,-2,-1,+0, -5.0e-01), -999.00000000, -0.00004, ""},
{/*185355*/DoodsonConstituent( 1,+3,+0,-2,+0,+0, +5.0e-01),  16.12989000, +0.00169, ""},
{/*185365*/DoodsonConstituent( 1,+3,+0,-2,+1,+0, +5.0e-01), -999.00000000, +0.00034, ""},
{/*185555*/DoodsonConstituent( 1,+3,+0,+0,+0,+0, +5.0e-01),  16.13910170, +0.01129, "Oo1"},
{/*185565*/DoodsonConstituent( 1,+3,+0,+0,+1,+0, +5.0e-01),  16.14131000, +0.00723, ""},
{/*185575*/DoodsonConstituent( 1,+3,+0,+0,+2,+0, +5.0e-01),  16.14352000, +0.00151, ""},
{/*185585*/DoodsonConstituent( 1,+3,+0,+0,+3,+0, +5.0e-01), -999.00000000, +0.00010, ""},
{/*186554*/DoodsonConstituent( 1,+3,+1,+0,+0,-1, -5.0e-01), -999.00000000, -0.00004, ""},
{/*191655*/DoodsonConstituent( 1,+4,-4,+1,+0,+0, +5.0e-01), -999.00000000, +0.00010, ""},
{/*192456*/DoodsonConstituent( 1,+4,-3,-1,+0,+1, +5.0e-01), -999.00000000, +0.00004, ""},
{/*193455*/DoodsonConstituent( 1,+4,-2,-1,+0,+0, +5.0e-01), -999.00000000, +0.00054, ""},
{/*193465*/DoodsonConstituent( 1,+4,-2,-1,+1,+0, +5.0e-01), -999.00000000, +0.00011, ""},
{/*193655*/DoodsonConstituent( 1,+4,-2,+1,+0,+0, +5.0e-01), -999.00000000, +0.00041, ""},
{/*193665*/DoodsonConstituent( 1,+4,-2,+1,+1,+0, +5.0e-01), -999.00000000, +0.00026, ""},
{/*193675*/DoodsonConstituent( 1,+4,-2,+1,+2,+0, +5.0e-01), -999.00000000, +0.00005, ""},
{/*195255*/DoodsonConstituent( 1,+4,+0,-3,+0,+0, +5.0e-01), -999.00000000, +0.00013, ""},
{/*195455*/DoodsonConstituent( 1,+4,+0,-1,+0,+0, +5.0e-01),  16.68347640, +0.00216, "ni1"},
{/*195465*/DoodsonConstituent( 1,+4,+0,-1,+1,+0, +5.0e-01),  16.68569000, +0.00138, ""},
{/*195475*/DoodsonConstituent( 1,+4,+0,-1,+2,+0, +5.0e-01), -999.00000000, +0.00029, ""},
//{/*209655*/DoodsonConstituent( 2,-5,+4,+1,+0,+0, -0.0e+00),  26.40793794, -999.00000, ""},
{/*215955*/DoodsonConstituent( 2,-4,+0,+4,+0,+0, +0.0e+00), -999.00000000, +0.00019, ""},
{/*217755*/DoodsonConstituent( 2,-4,+2,+2,+0,+0, +0.0e+00),  26.87945903, +0.00078, ""},
{/*218754*/DoodsonConstituent( 2,-4,+3,+2,+0,-1, +0.0e+00), -999.00000000, +0.00006, ""},
{/*219555*/DoodsonConstituent( 2,-4,+4,+0,+0,+0, +0.0e+00),  26.95231264, +0.00048, ""},
//{/*219755*/DoodsonConstituent( 2,-4,+4,+2,+0,+0, -0.0e+00),  26.96159631, -999.00000, ""},
//{/*225655*/DoodsonConstituent( 2,-3,+0,+1,+0,+0, +1.0e+00),  27.34169645, -999.00000, ""},
{/*225845*/DoodsonConstituent( 2,-3,+0,+3,-1,+0, +1.0e+00), -999.00000000, -0.00007, ""},
{/*225855*/DoodsonConstituent( 2,-3,+0,+3,+0,+0, +0.0e+00), -999.00000000, +0.00180, "3N2"},
{/*226656*/DoodsonConstituent( 2,-3,+1,+1,+0,+1, +1.0e+00), -999.00000000, -0.00009, ""},
{/*226854*/DoodsonConstituent( 2,-3,+1,+3,+0,-1, +0.0e+00), -999.00000000, +0.00004, ""},
{/*227645*/DoodsonConstituent( 2,-3,+2,+1,-1,+0, +1.0e+00), -999.00000000, -0.00017, ""},
{/*227655*/DoodsonConstituent( 2,-3,+2,+1,+0,+0, +0.0e+00),  27.42383373, +0.00467, "eps2"},
{/*228654*/DoodsonConstituent( 2,-3,+3,+1,+0,-1, +0.0e+00), -999.00000000, +0.00036, ""},
{/*229445*/DoodsonConstituent( 2,-3,+4,-1,-1,+0, +1.0e+00), -999.00000000, -0.00003, ""},
{/*229455*/DoodsonConstituent( 2,-3,+4,-1,+0,+0, +0.0e+00),  27.49668734, +0.00090, ""},
//{/*229655*/DoodsonConstituent( 2,-3,+4,+1,+0,+0, -0.0e+00),  27.50597101, -999.00000, ""},
//{/*233555*/DoodsonConstituent( 2,-2,-2,+0,+0,+0, -0.0e+00),  27.80393387, -999.00000, ""},
{/*233955*/DoodsonConstituent( 2,-2,-2,+4,+0,+0, +1.0e+00), -999.00000000, -0.00006, ""},
{/*234656*/DoodsonConstituent( 2,-2,-1,+1,+0,+1, +0.0e+00), -999.00000000, +0.00042, ""},
{/*234756*/DoodsonConstituent( 2,-2,-1,+2,+0,+1, +1.0e+00), -999.00000000, -0.00022, ""},
{/*235535*/DoodsonConstituent( 2,-2,+0,+0,-2,+0, +1.0e+00), -999.00000000, -0.00010, ""},
//{/*235555*/DoodsonConstituent( 2,-2,+0,+0,+0,+0, -0.0e+00),  27.88607115, -999.00000, ""},
{/*235745*/DoodsonConstituent( 2,-2,+0,+2,-1,+0, +1.0e+00), -999.00000000, -0.00060, ""},
{/*235755*/DoodsonConstituent( 2,-2,+0,+2,+0,+0, +0.0e+00),  27.89535482, +0.01601, "2N2"},
{/*236556*/DoodsonConstituent( 2,-2,+1,+0,+0,+1, +1.0e+00), -999.00000000, -0.00027, ""},
{/*236655*/DoodsonConstituent( 2,-2,+1,+1,+0,+0, +1.0e+00), -999.00000000, -0.00017, ""},
{/*236754*/DoodsonConstituent( 2,-2,+1,+2,+0,-1, +0.0e+00), -999.00000000, +0.00025, ""},
{/*237545*/DoodsonConstituent( 2,-2,+2,+0,-1,+0, +1.0e+00), -999.00000000, -0.00072, ""},
{/*237555*/DoodsonConstituent( 2,-2,+2,+0,+0,+0, +0.0e+00),  27.96820843, +0.01932, "mu2"},
{/*238455*/DoodsonConstituent( 2,-2,+3,-1,+0,+0, +1.0e+00), -999.00000000, -0.00004, ""},
{/*238544*/DoodsonConstituent( 2,-2,+3,+0,-1,-1, +1.0e+00), -999.00000000, -0.00005, ""},
{/*238554*/DoodsonConstituent( 2,-2,+3,+0,+0,-1, +0.0e+00), -999.00000000, +0.00130, ""},
{/*239355*/DoodsonConstituent( 2,-2,+4,-2,+0,+0, +0.0e+00), -999.00000000, +0.00059, ""},
{/*239553*/DoodsonConstituent( 2,-2,+4,+0,+0,-2, +0.0e+00), -999.00000000, +0.00005, ""},
{/*243635*/DoodsonConstituent( 2,-1,-2,+1,-2,+0, +1.0e+00), -999.00000000, -0.00010, ""},
//{/*243655*/DoodsonConstituent( 2,-1,-2,+1,+0,+0, -0.0e+00),  28.35759224, -999.00000, ""},
{/*243855*/DoodsonConstituent( 2,-1,-2,+3,+0,+0, +1.0e+00), -999.00000000, -0.00039, ""},
{/*244646*/DoodsonConstituent( 2,-1,-1,+1,-1,+1, +0.0e+00), -999.00000000, +0.00003, ""},
{/*244656*/DoodsonConstituent( 2,-1,-1,+1,+0,+1, +1.0e+00), -999.00000000, -0.00102, ""},
{/*245435*/DoodsonConstituent( 2,-1,+0,-1,-2,+0, +1.0e+00), -999.00000000, -0.00047, ""},
{/*245555*/DoodsonConstituent( 2,-1,+0,+0,+0,+0, +1.0e+00), -999.00000000, -0.00389, ""},
{/*245556*/DoodsonConstituent( 2,-1,+0,+0,+0,+1, +0.0e+00), -999.00000000, +0.00010, ""},
{/*245635*/DoodsonConstituent( 2,-1,+0,+1,-2,+0, +0.0e+00), -999.00000000, +0.00007, ""},
{/*245645*/DoodsonConstituent( 2,-1,+0,+1,-1,+0, +1.0e+00), -999.00000000, -0.00451, ""},
{/*245655*/DoodsonConstituent( 2,-1,+0,+1,+0,+0, +0.0e+00),  28.43972952, +0.12099, "N2"},
{/*246456*/DoodsonConstituent( 2,-1,+1,-1,+0,+1, +1.0e+00), -999.00000000, -0.00022, ""},
{/*246555*/DoodsonConstituent( 2,-1,+1,+0,+0,+0, +1.0e+00), -999.00000000, -0.00065, ""},
{/*246644*/DoodsonConstituent( 2,-1,+1,+1,-1,-1, +1.0e+00), -999.00000000, -0.00004, ""},
{/*246654*/DoodsonConstituent( 2,-1,+1,+1,+0,-1, +0.0e+00), -999.00000000, +0.00113, ""},
{/*247445*/DoodsonConstituent( 2,-1,+2,-1,-1,+0, +1.0e+00), -999.00000000, -0.00086, ""},
{/*247455*/DoodsonConstituent( 2,-1,+2,-1,+0,+0, +0.0e+00),  28.51258312, +0.02298, "ni2"},
{/*247655*/DoodsonConstituent( 2,-1,+2,+1,+0,+0, +0.0e+00), -999.00000000, +0.00010, ""},
{/*247665*/DoodsonConstituent( 2,-1,+2,+1,+1,+0, +1.0e+00), -999.00000000, -0.00008, ""},
{/*248444*/DoodsonConstituent( 2,-1,+3,-1,-1,-1, +1.0e+00), -999.00000000, -0.00004, ""},
{/*248454*/DoodsonConstituent( 2,-1,+3,-1,+0,-1, +0.0e+00), -999.00000000, +0.00106, ""},
{/*252756*/DoodsonConstituent( 2,+0,-3,+2,+0,+1, +1.0e+00), -999.00000000, -0.00008, ""},
{/*253535*/DoodsonConstituent( 2,+0,-2,+0,-2,+0, +1.0e+00), -999.00000000, -0.00028, ""},
//{/*253555*/DoodsonConstituent( 2,+0,-2,+0,+0,+0, -0.0e+00),  28.90196694, -999.00000, ""},
{/*253745*/DoodsonConstituent( 2,+0,-2,+2,-1,+0, +0.0e+00), -999.00000000, +0.00007, ""},
{/*253755*/DoodsonConstituent( 2,+0,-2,+2,+0,+0, +1.0e+00), -999.00000000, -0.00190, "gamma2"},
{/*254546*/DoodsonConstituent( 2,+0,-1,+0,-1,+1, +0.0e+00), -999.00000000, +0.00005, ""},
//{/*254555*/DoodsonConstituent( 2,+0,-1,+0,+0,+0, +5.0e-01),  28.94303557, -999.00000, ""},
{/*254556*/DoodsonConstituent( 2,+0,-1,+0,+0,+1, +1.0e+00), -999.00000000, -0.00218, "alpha2"},
{/*254655*/DoodsonConstituent( 2,+0,-1,+1,+0,+0, +0.0e+00), -999.00000000, +0.00009, ""},
{/*255535*/DoodsonConstituent( 2,+0,+0,+0,-2,+0, +0.0e+00), -999.00000000, +0.00033, ""},
{/*255545*/DoodsonConstituent( 2,+0,+0,+0,-1,+0, +1.0e+00), -999.00000000, -0.02358, ""},
{/*255555*/DoodsonConstituent( 2,+0,+0,+0,+0,+0, +0.0e+00),  28.98410421, +0.63192, "M2"},
{/*255755*/DoodsonConstituent( 2,+0,+0,+2,+0,+0, +0.0e+00), -999.00000000, +0.00037, ""},
{/*255765*/DoodsonConstituent( 2,+0,+0,+2,+1,+0, +0.0e+00), -999.00000000, +0.00013, ""},
{/*256544*/DoodsonConstituent( 2,+0,+1,+0,-1,-1, +1.0e+00), -999.00000000, -0.00004, ""},
{/*256554*/DoodsonConstituent( 2,+0,+1,+0,+0,-1, +0.0e+00), -999.00000000, +0.00192, "beta2"},
//{/*256555*/DoodsonConstituent( 2,+0,+1,+0,+0,+0, -5.0e-01),  29.02517285, -999.00000, ""},
{/*257355*/DoodsonConstituent( 2,+0,+2,-2,+0,+0, +1.0e+00), -999.00000000, -0.00036, ""},
{/*257555*/DoodsonConstituent( 2,+0,+2,+0,+0,+0, +0.0e+00),  29.06624149, +0.00072, "dlt2"},
{/*257565*/DoodsonConstituent( 2,+0,+2,+0,+1,+0, +1.0e+00), -999.00000000, -0.00036, ""},
{/*257575*/DoodsonConstituent( 2,+0,+2,+0,+2,+0, +0.0e+00), -999.00000000, +0.00012, ""},
{/*258554*/DoodsonConstituent( 2,+0,+3,+0,+0,-1, +0.0e+00), -999.00000000, +0.00005, ""},
{/*262656*/DoodsonConstituent( 2,+1,-3,+1,+0,+1, +1.0e+00), -999.00000000, -0.00022, ""},
{/*263645*/DoodsonConstituent( 2,+1,-2,+1,-1,+0, +0.0e+00), -999.00000000, +0.00021, ""},
{/*263655*/DoodsonConstituent( 2,+1,-2,+1,+0,+0, +1.0e+00),  29.45562530, -0.00466, "lambda2"},
{/*264456*/DoodsonConstituent( 2,+1,-1,-1,+0,+1, +1.0e+00), -999.00000000, -0.00007, ""},
{/*264555*/DoodsonConstituent( 2,+1,-1,+0,+0,+0, +0.0e+00), -999.00000000, +0.00011, ""},
{/*265445*/DoodsonConstituent( 2,+1,+0,-1,-1,+0, +0.0e+00), -999.00000000, +0.00066, ""},
{/*265455*/DoodsonConstituent( 2,+1,+0,-1,+0,+0, +1.0e+00),  29.52847891, -0.01786, "L2"},
{/*265555*/DoodsonConstituent( 2,+1,+0,+0,+0,+0, +0.0e+00), -999.00000000, +0.00359, ""},
{/*265645*/DoodsonConstituent( 2,+1,+0,+1,-1,+0, +1.0e+00), -999.00000000, -0.00008, ""},
{/*265655*/DoodsonConstituent( 2,+1,+0,+1,+0,+0, +0.0e+00),  29.53776258, +0.00447, ""},
{/*265665*/DoodsonConstituent( 2,+1,+0,+1,+1,+0, +0.0e+00), -999.00000000, +0.00197, ""},
{/*265675*/DoodsonConstituent( 2,+1,+0,+1,+2,+0, +0.0e+00), -999.00000000, +0.00028, ""},
{/*267455*/DoodsonConstituent( 2,+1,+2,-1,+0,+0, +0.0e+00), -999.00000000, +0.00086, ""},
{/*267465*/DoodsonConstituent( 2,+1,+2,-1,+1,+0, +0.0e+00), -999.00000000, +0.00041, ""},
{/*267475*/DoodsonConstituent( 2,+1,+2,-1,+2,+0, +0.0e+00), -999.00000000, +0.00005, ""},
{/*271557*/DoodsonConstituent( 2,+2,-4,+0,+0,+2, +0.0e+00), -999.00000000, +0.00070, "2t2"},
{/*272556*/DoodsonConstituent( 2,+2,-3,+0,+0,+1, +0.0e+00),  29.95893332, +0.01720, "T2"},
{/*273545*/DoodsonConstituent( 2,+2,-2,+0,-1,+0, +0.0e+00), -999.00000000, +0.00066, ""},
{/*273555*/DoodsonConstituent( 2,+2,-2,+0,+0,+0, +0.0e+00),  30.00000000, +0.29400, "S2"},
{/*273755*/DoodsonConstituent( 2,+2,-2,+2,+0,+0, +0.0e+00), -999.00000000, +0.00004, ""},
{/*274554*/DoodsonConstituent( 2,+2,-1,+0,+0,-1, +1.0e+00),  30.04106668, -0.00246, "R2"},
{/*274556*/DoodsonConstituent( 2,+2,-1,+0,+0,+1, +0.0e+00), -999.00000000, +0.00062, ""},
{/*274566*/DoodsonConstituent( 2,+2,-1,+0,+1,+1, +1.0e+00), -999.00000000, -0.00004, ""},
{/*275545*/DoodsonConstituent( 2,+2,+0,+0,-1,+0, +1.0e+00), -999.00000000, -0.00102, ""},
{/*275555*/DoodsonConstituent( 2,+2,+0,+0,+0,+0, +0.0e+00),  30.08213728, +0.07996, "K2"},
{/*275565*/DoodsonConstituent( 2,+2,+0,+0,+1,+0, +0.0e+00), -999.00000000, +0.02383, "K2+"},
{/*275575*/DoodsonConstituent( 2,+2,+0,+0,+2,+0, +0.0e+00), -999.00000000, +0.00259, "K2++"},
{/*276554*/DoodsonConstituent( 2,+2,+1,+0,+0,-1, +0.0e+00), -999.00000000, +0.00063, ""},
{/*277355*/DoodsonConstituent( 2,+2,+2,-2,+0,+0, +0.0e+00), -999.00000000, +0.00004, ""},
{/*277555*/DoodsonConstituent( 2,+2,+2,+0,+0,+0, +0.0e+00), -999.00000000, +0.00053, ""},
{/*282656*/DoodsonConstituent( 2,+3,-3,+1,+0,+1, +0.0e+00), -999.00000000, +0.00004, ""},
{/*283445*/DoodsonConstituent( 2,+3,-2,-1,-1,+0, +0.0e+00), -999.00000000, +0.00006, ""},
{/*283455*/DoodsonConstituent( 2,+3,-2,-1,+0,+0, +0.0e+00),  30.54437470, +0.00004, ""},
{/*283655*/DoodsonConstituent( 2,+3,-2,+1,+0,+0, +0.0e+00), -999.00000000, +0.00086, "ksi2"},
{/*283665*/DoodsonConstituent( 2,+3,-2,+1,+1,+0, +0.0e+00), -999.00000000, +0.00037, ""},
{/*283675*/DoodsonConstituent( 2,+3,-2,+1,+2,+0, +0.0e+00), -999.00000000, +0.00004, ""},
{/*285445*/DoodsonConstituent( 2,+3,+0,-1,-1,+0, +1.0e+00), -999.00000000, -0.00009, ""},
{/*285455*/DoodsonConstituent( 2,+3,+0,-1,+0,+0, +0.0e+00),  30.62651197, +0.00447, "ita2"},
{/*285465*/DoodsonConstituent( 2,+3,+0,-1,+1,+0, +0.0e+00), -999.00000000, +0.00195, ""},
{/*285475*/DoodsonConstituent( 2,+3,+0,-1,+2,+0, +0.0e+00), -999.00000000, +0.00022, ""},
{/*285655*/DoodsonConstituent( 2,+3,+0,+1,+0,+0, +1.0e+00), -999.00000000, -0.00003, ""},
//{/*291555*/DoodsonConstituent( 2,+4,-4,+0,+0,+0, -0.0e+00),  31.01589579, -999.00000, ""},
{/*292556*/DoodsonConstituent( 2,+4,-3,+0,+0,+1, +0.0e+00), -999.00000000, +0.00005, ""},
{/*293555*/DoodsonConstituent( 2,+4,-2,+0,+0,+0, +0.0e+00),  31.09803306, +0.00074, ""},
{/*293565*/DoodsonConstituent( 2,+4,-2,+0,+1,+0, +0.0e+00), -999.00000000, +0.00032, ""},
{/*293575*/DoodsonConstituent( 2,+4,-2,+0,+2,+0, +0.0e+00), -999.00000000, +0.00003, ""},
{/*295355*/DoodsonConstituent( 2,+4,+0,-2,+0,+0, +0.0e+00), -999.00000000, +0.00037, ""},
{/*295365*/DoodsonConstituent( 2,+4,+0,-2,+1,+0, +0.0e+00), -999.00000000, +0.00016, ""},
{/*295555*/DoodsonConstituent( 2,+4,+0,+0,+0,+0, +0.0e+00), -999.00000000, +0.00117, ""},
{/*295565*/DoodsonConstituent( 2,+4,+0,+0,+1,+0, +0.0e+00), -999.00000000, +0.00101, ""},
{/*295575*/DoodsonConstituent( 2,+4,+0,+0,+2,+0, +0.0e+00), -999.00000000, +0.00033, ""},
{/*295585*/DoodsonConstituent( 2,+4,+0,+0,+3,+0, +0.0e+00), -999.00000000, +0.00005, ""},
//{/*301455*/DoodsonConstituent( 3,-5,-4,-1,+0,+0, -0.0e+00),  31.56027048, -999.00000, ""},
//{/*309555*/DoodsonConstituent( 3,-5,+4,+0,+0,+0, -0.0e+00),  32.03179157, -999.00000, ""},
//{/*335655*/DoodsonConstituent( 3,-2,+0,+1,+0,+0, +5.0e-01),  42.38276509, -999.00000, ""},
//{/*345555*/DoodsonConstituent( 3,-1,+0,+0,+0,+0, +5.0e-01),  42.92713979, -999.00000, ""},
//{/*355555*/DoodsonConstituent( 3,+0,+0,+0,+0,+0, +1.0e+00),  43.47615632, -999.00000, "m3"},
//{/*363555*/DoodsonConstituent( 3,+1,-2,+0,+0,+0, +5.0e-01),  43.94303558, -999.00000, ""},
//{/*364555*/DoodsonConstituent( 3,+1,-1,+0,+0,+0, -0.0e+00),  43.98410421, -999.00000, ""},
//{/*365555*/DoodsonConstituent( 3,+1,+0,+0,+0,+0, -5.0e-01),  44.02517285, -999.00000, ""},
//{/*381555*/DoodsonConstituent( 3,+3,-4,+0,+0,+0, +1.0e+00),  44.95893136, -999.00000, "T3"},
//{/*382555*/DoodsonConstituent( 3,+3,-3,+0,+0,+0, +1.0e+00),  45.00000000, -999.00000, "s3"},
//{/*383555*/DoodsonConstituent( 3,+3,-2,+0,+0,+0, +1.0e+00),  45.04106864, -999.00000, "r3"},
//{/*385555*/DoodsonConstituent( 3,+3,+0,+0,+0,+0, -1.5e+00),  45.12320592, -999.00000, ""},
//{/*427655*/DoodsonConstituent( 4,-3,+2,+1,+0,+0, -0.0e+00),  56.40793794, -999.00000, ""},
//{/*435755*/DoodsonConstituent( 4,-2,+0,+2,+0,+0, +0.0e+00),  56.87945903, -999.00000, "n4"},
//{/*437555*/DoodsonConstituent( 4,-2,+2,+0,+0,+0, -0.0e+00),  56.95231264, -999.00000, ""},
//{/*445655*/DoodsonConstituent( 4,-1,+0,+1,+0,+0, +0.0e+00),  57.42383373, -999.00000, "mn4"},
//{/*447455*/DoodsonConstituent( 4,-1,+2,-1,+0,+0, -0.0e+00),  57.49668734, -999.00000, ""},
//{/*455555*/DoodsonConstituent( 4,+0,+0,+0,+0,+0, +0.0e+00),  57.96820843, -999.00000, "m4"},
//{/*463655*/DoodsonConstituent( 4,+1,-2,+1,+0,+0, -0.0e+00),  58.43972952, -999.00000, ""},
//{/*465455*/DoodsonConstituent( 4,+1,+0,-1,+0,+0, +1.0e+00),  58.51258312, -999.00000, ""},
//{/*465655*/DoodsonConstituent( 4,+1,+0,+1,+0,+0, -0.0e+00),  58.52186679, -999.00000, ""},
//{/*473555*/DoodsonConstituent( 4,+2,-2,+0,+0,+0, +0.0e+00),  58.98410421, -999.00000, "ms4"},
//{/*475555*/DoodsonConstituent( 4,+2,+0,+0,+0,+0, -0.0e+00),  59.06624149, -999.00000, ""},
//{/*483455*/DoodsonConstituent( 4,+3,-2,-1,+0,+0, -0.0e+00),  59.52847891, -999.00000, ""},
//{/*491555*/DoodsonConstituent( 4,+4,-4,+0,+0,+0, +0.0e+00),  60.00000000, -999.00000, "s4"},
//{/*493555*/DoodsonConstituent( 4,+4,-2,+0,+0,+0, -0.0e+00),  60.08213728, -999.00000, ""},
//{/*625655*/DoodsonConstituent( 6,-3,+0,+1,+0,+0, -0.0e+00),  85.30990488, -999.00000, ""},
//{/*627655*/DoodsonConstituent( 6,-3,+2,+1,+0,+0, -0.0e+00),  85.39204216, -999.00000, ""},
//{/*635555*/DoodsonConstituent( 6,-2,+0,+0,+0,+0, -0.0e+00),  85.85427958, -999.00000, ""},
//{/*635755*/DoodsonConstituent( 6,-2,+0,+2,+0,+0, -1.0e+00),  85.86356325, -999.00000, ""},
//{/*637555*/DoodsonConstituent( 6,-2,+2,+0,+0,+0, -0.0e+00),  85.93641686, -999.00000, ""},
//{/*645655*/DoodsonConstituent( 6,-1,+0,+1,+0,+0, -0.0e+00),  86.40793794, -999.00000, ""},
//{/*647455*/DoodsonConstituent( 6,-1,+2,-1,+0,+0, -0.0e+00),  86.48079155, -999.00000, ""},
//{/*650435*/DoodsonConstituent( 6,+0,-5,-1,-2,+0, -1.0e+03), -999.00000000, -0.00003, ""},
//{/*653555*/DoodsonConstituent( 6,+0,-2,+0,+0,+0, -0.0e+00),  86.87017536, -999.00000, ""},
//{/*655555*/DoodsonConstituent( 6,+0,+0,+0,+0,+0, +0.0e+00),  86.95231264, -999.00000, "m6"},
//{/*657555*/DoodsonConstituent( 6,+0,+2,+0,+0,+0, -0.0e+00),  87.03444992, -999.00000, ""},
//{/*663655*/DoodsonConstituent( 6,+1,-2,+1,+0,+0, -0.0e+00),  87.42383373, -999.00000, ""},
//{/*665455*/DoodsonConstituent( 6,+1,+0,-1,+0,+0, +1.0e+00),  87.49668734, -999.00000, ""},
//{/*665655*/DoodsonConstituent( 6,+1,+0,+1,+0,+0, -0.0e+00),  87.50597101, -999.00000, ""},
//{/*673555*/DoodsonConstituent( 6,+2,-2,+0,+0,+0, -0.0e+00),  87.96820843, -999.00000, ""},
//{/*675355*/DoodsonConstituent( 6,+2,+0,-2,+0,+0, +1.0e+00),  88.04106204, -999.00000, ""},
//{/*675555*/DoodsonConstituent( 6,+2,+0,+0,+0,+0, -0.0e+00),  88.05034571, -999.00000, ""},
//{/*683455*/DoodsonConstituent( 6,+3,-2,-1,+0,+0, -0.0e+00),  88.51258313, -999.00000, ""},
//{/*685455*/DoodsonConstituent( 6,+3,+0,-1,+0,+0, -0.0e+00),  88.59472040, -999.00000, ""},
//{/*691555*/DoodsonConstituent( 6,+4,-4,+0,+0,+0, -0.0e+00),  88.98410421, -999.00000, ""},
//{/*693555*/DoodsonConstituent( 6,+4,-2,+0,+0,+0, -0.0e+00),  89.06624149, -999.00000, ""},
//{/*800656*/DoodsonConstituent( 8,-5,-5,+1,+0,+1, -1.0e+03), -999.00000000, -0.00002, ""},
//{/*835755*/DoodsonConstituent( 8,-2,+0,+2,+0,+0, -0.0e+00), 114.84766746, -999.00000, ""},
//{/*845655*/DoodsonConstituent( 8,-1,+0,+1,+0,+0, -0.0e+00), 115.39204216, -999.00000, ""},
//{/*850455*/DoodsonConstituent( 8,+0,-5,-1,+0,+0, -1.0e+03), -999.00000000, -0.01276, ""},
//{/*855555*/DoodsonConstituent( 8,+0,+0,+0,+0,+0, +0.0e+00), 115.93641686, -999.00000, "m8"},
//{/*860464*/DoodsonConstituent( 8,+1,-5,-1,+1,-1, -1.0e+03), -999.00000000, +0.00004, ""},
//{/*863655*/DoodsonConstituent( 8,+1,-2,+1,+0,+0, -0.0e+00), 116.40793794, -999.00000, ""},
//{/*865455*/DoodsonConstituent( 8,+1,+0,-1,+0,+0, +1.0e+00), 116.48079155, -999.00000, ""},
//{/*873555*/DoodsonConstituent( 8,+2,-2,+0,+0,+0, -0.0e+00), 116.95231264, -999.00000, ""},
//{/*875555*/DoodsonConstituent( 8,+2,+0,+0,+0,+0, -0.0e+00), 117.03444992, -999.00000, ""},
//{/*891555*/DoodsonConstituent( 8,+4,-4,+0,+0,+0, -0.0e+00), 117.96820843, -999.00000, ""},
//{/*893555*/DoodsonConstituent( 8,+4,-2,+0,+0,+0, -0.0e+00), 118.05034571, -999.00000, ""},
//{/*895555*/DoodsonConstituent( 8,+4,+0,+0,+0,+0, -0.0e+00), 118.13248298, -999.00000, ""},
//{/*950355*/DoodsonConstituent( 9,+0,-5,-2,+0,+0, -1.0e+03), -999.00000000, -0.00169, ""},
    }}; /* TidalConstituentsArray */
} /* namespace dso */
#endif
