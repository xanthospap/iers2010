/** @file
 * Define a DoodsonConstituent class to represent tidal constituents via
 * Doodson numbers.
 *
 * References:
 * [1]
 * https://ivscc.gsfc.nasa.gov/hfeop_wg/memos/memo-conventions_Ray_2017Dec10.pdf
 *
 * [2] Balidakis, K., Sulzbach, R., Shihora, L., Dahle, C., Dill, R., &
 * Dobslaw, H. (2022). Atmospheric contributions to global ocean tides for
 * satellite gravimetry. Journal of Advances in Modeling Earth Systems, 14,
 * e2022MS003193. https://doi.org/10.1029/2022MS003193
 *
 * [3] Henryk Dobslaw, Inga Bergmann-Wolf, Robert Dill, Lea Poropat, Frank
 * Flechtner, GRACE 327-750 Gravity Recovery and Climate Experiment Product
 * Description Document for AOD1B Release 06 (Rev. 6.1, October 19, 2017),
 * GFZ German Research Centre for Geosciences Department 1: Geodesy
 *
 * [4] Daniel Rieser, Torsten Mayer-G¨urr, Roman Savcenko, Wolfgang Bosch, 
 * Johann W¨unsch, Christoph Dahle and Frank Flechtner, The ocean tide model 
 * EOT11a in spherical harmonics representation, Technical Note, July 2012,
 * https://www.tugraz.at/fileadmin/user_upload/Institute/IFG/satgeo/pdf/TN_EOT11a.pdf
 *
 * [5] Torsten Mayer-Guerr, 00README_simulation.txt, available from the data 
 * download area for the COST-G benchmark test, see [6]
 *
 * [6] Lasser, M., Meyer, U., Jäggi, A., Mayer-Gürr, T., Kvas, A., Neumayer,
 * K. H., Dahle, C., Flechtner, F., Lemoine, J.-M., Koch, I., Weigelt, M.,
 * and Flury, J.: Benchmark data for verifying background model
 * implementations in orbit and gravity field determination software,
 * Adv. Geosci., 55, 1–11, https://doi.org/10.5194/adgeo-55-1-2020, 2020.
 */

#ifndef __DSO_DOODSON_NUMBER_DEFINES_HPP__
#define __DSO_DOODSON_NUMBER_DEFINES_HPP__

#include "datetime/calendar.hpp"
#include "geodesy/units.hpp"
#include <array>
#include <cassert>
#include <cctype>
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

  /** Convert a character (part of a Doodson number string) to int.
   *
   * The conversion follows the GROOPS convention, which is:
   *   [n] -> -13 [0] -> 0 [a] -> 10
   *   [o] -> -12 [1] -> 1 [b] -> 11
   *   [p] -> -11 [2] -> 2 [c] -> 12
   *   [q] -> -10 [3] -> 3 [d] -> 13
   *   [r] -> -9  [4] -> 4 [e] -> 14
   *   [s] -> -8  [5] -> 5 [f] -> 15
   *   [t] -> -7  [6] -> 6 [g] -> 16
   *   [u] -> -6  [7] -> 7 [h] -> 17
   *   [v] -> -5  [8] -> 8 [i] -> 18
   *   [x] -> -3  [9] -> 9 [j] -> 19
   *   [y] -> -2           [k] -> 20
   *   [z] -> -1           [l] -> 21
   *   [a] -> 10           [m] -> 22
   *
   *  The function is case-insesitive, i.e. upper- and lower-case characters
   *  are accepted.
   *  Note that this conversion ignores the '+5' convension usually followed
   *  when resolving Doodson numbers. If needed, subtract 5 from the resolved
   *  int.
   */
  inline static constexpr int char2int(char c) {
    if (std::isdigit(c))
      return (c - '0');
    c = std::tolower(c);
    if (('a' <= c) && (c <= 'm'))
      return (c - 'a') + 10;
    if (('n' <= c) && (c <= 'z'))
      return (c - 'n') - 13;
    throw std::runtime_error(
        "[ERROR] Failed converting char to doodson integer\n");
  }

  /** Convert an int (part of a Doodson number) to character.
   *
   * The conversion follows the GROOPS convention, which is:
   *   [n] -> -13 [0] -> 0 [a] -> 10
   *   [o] -> -12 [1] -> 1 [b] -> 11
   *   [p] -> -11 [2] -> 2 [c] -> 12
   *   [q] -> -10 [3] -> 3 [d] -> 13
   *   [r] -> -9  [4] -> 4 [e] -> 14
   *   [s] -> -8  [5] -> 5 [f] -> 15
   *   [t] -> -7  [6] -> 6 [g] -> 16
   *   [u] -> -6  [7] -> 7 [h] -> 17
   *   [v] -> -5  [8] -> 8 [i] -> 18
   *   [x] -> -3  [9] -> 9 [j] -> 19
   *   [y] -> -2           [k] -> 20
   *   [z] -> -1           [l] -> 21
   *   [a] -> 10           [m] -> 22
   */
  inline static constexpr char int2char(int d) noexcept {
    if (d < 0)
      return ('n' + (d + 13));
    if (d > 9)
      return ('a' + (d - 10));
    return ('0' + d);
  }

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

  /** @brief Resolve a string of type: '272.556' to a Doodson number.
   *
   * This is how waves are usually denoted in the IERS 2010 standards, using
   * the Doodson convention.
   *
   * The input string does not have to be null-terminated. Only the first
   * 7 characters will be considered.
   *
   * Note that IERS 2010 Doodson strings comply with the '+5' convnention,
   * i.e. the actual multipliers are the integers given minus 5 (except for
   * the first one). Here, this convention is followed, hence e.g. the Doodson
   * number 123.456, will resolve to integers {1,-3,-2, -1, 0, 1}.
   */
  static DoodsonConstituent from_chars(const char *str);

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
                               double pifactor = 0) noexcept {
    iar[0] = a0;
    iar[1] = a1;
    iar[2] = a2;
    iar[3] = a3;
    iar[4] = a4;
    iar[5] = a5;
    pifac = pifactor;
  }

  /* @brief Transform to a Doodson-number string as: "xxx.xxx"
   *
   * For integers that are <0 or >9, we follow the GROOPS convention, i.e.:
   *   [n] -> -13 [0] -> 0 [a] -> 10
   *   [o] -> -12 [1] -> 1 [b] -> 11
   *   [p] -> -11 [2] -> 2 [c] -> 12
   *   [q] -> -10 [3] -> 3 [d] -> 13
   *   [r] -> -9  [4] -> 4 [e] -> 14
   *   [s] -> -8  [5] -> 5 [f] -> 15
   *   [t] -> -7  [6] -> 6 [g] -> 16
   *   [u] -> -6  [7] -> 7 [h] -> 17
   *   [v] -> -5  [8] -> 8 [i] -> 18
   *   [x] -> -3  [9] -> 9 [j] -> 19
   *   [y] -> -2           [k] -> 20
   *   [z] -> -1           [l] -> 21
   *   [a] -> 10           [m] -> 22
   *
   * @param[out] buf A char buffer of size at least 7 where the instance
   *            will be written to.
   * @param[in] use_5s_convention Add +5 to all integers after the first
   *            place (i.e. write the string as: τ,s+5,h+5,p+5,N'+5,ps+5)
   * @return A pointer to the start of the string, i.e. to buf
   */
  char *str(char *buf, bool use_5s_convention = true) const noexcept;

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

  double pifactor() const noexcept { return pifac; }

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
  double argument(const double *const doodson_arguments) const noexcept {
    const double *__restrict__ f = doodson_arguments;
    return dso::anp(f[0] * iar[0] + f[1] * iar[1] + f[2] * iar[2] +
                    f[3] * iar[3] + f[4] * iar[4] + f[5] * iar[5]);
  }

  const int *int_array() const noexcept { return iar; }

}; /* DoodsonConstituent */

/** @brief Resolve a string of type: '272.556' to a Doodson number/wave.
 *
 * This is how waves are usually denoted in the IERS 2010 standards, using
 * the Doodson convention.
 *
 * The input string does not have to be null-terminated. Only the first
 * 7 characters will be considered, AFTER ommiting any leading whitespace
 * characters. Hence, the string can start with any number of whitespaces
 * wchich will be skipped.
 *
 * Note that IERS 2010 Dooson strings comply with the  '+5' convnention,
 * i.e. the actual multipliers are the integers given minus 5 (except for
 * the first one).
 */
DoodsonConstituent resolve_iers10_doodson_string(const char *);

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
 *   * [5] pl : Longitude of Sun's mean perigee [rad]; often also noted as p_s
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
/** @class
 * This class is the same as TidalWave with the only difference of a const
 * char member. This allows for constexpr construction, and is therefor
 * only used for constructing compile-time TidalWave's.
 * An instance of this class can be used to construct an instance of type
 * TidalWave.
 */
struct TidalConstituentArrayEntry {
  DoodsonConstituent _d;
  double _per; /* period in [deg/hour] */
  double _hf;  /* from IERS 2010, Table 6.7 and [2] */
  char _n[16];  /* name */
  const DoodsonConstituent &doodson() const noexcept { return _d; }
  DoodsonConstituent &doodson() noexcept { return _d; }
  double period() const noexcept { return _per; }
  double &period() noexcept { return _per; }
  double height() const noexcept { return _hf; }
  double &height() noexcept { return _hf; }
  const char *name() const noexcept { return _n; }
}; /* TidalConstituentsArrayEntry */

/* Atmospheric tidal waves as described in [2].
 * Note that we have no information on the Dooson-Warburg correction here.
 */
constexpr static std::array<TidalConstituentArrayEntry, 18>
    AtmosphericTidalHarmonics = {
        {{/*162.556*/ DoodsonConstituent(1, 1, -3, 0, 0, 1, -1.0e+00),
          14.917866090000, 0.00, "π1"},
         {/*163.555*/ DoodsonConstituent(1, 1, -2, 0, 0, 0, -1.0e+00),
          14.958932770000, 0.00, "P1"},
         // as in [2]
         //{/*164.556*/ DoodsonConstituent(1, 1, -1, 0, 0, 1, +1.0e+00),
         // 15.000001410000, 0.00, "S1"},
         /* as in [3] */
         {/*164.556*/ DoodsonConstituent(1, 1, -1, 0, 0, 0, +2.0e+00),
          15.000001410000, 0.00, "S1"},
         {/*165.555*/ DoodsonConstituent(1, 1, 0, 0, 0, 0, +1.0e+00),
          15.041070050000, 0.00, "K1"},
         {/*166.554*/ DoodsonConstituent(1, 1, 1, 0, 0, -1, +1.0e+00),
          15.082136730000, 0.00, "ψ1"},
         {/*255.555*/ DoodsonConstituent(2, 0, 0, 0, 0, 0, +0.0e+00),
          28.984107050000, 0.00, "M2"},
         {/*272.556*/ DoodsonConstituent(2, 2, -3, 0, 0, 1, +0.0e+00),
          29.958936140000, 0.00, "T2"},
         {/*273.555*/ DoodsonConstituent(2, 2, -2, 0, 0, 0, -1.0e+00),
          30.000002820000, 0.00, "S2"},
         {/*274.554*/ DoodsonConstituent(2, 2, -1, 0, 0, -1, +2.0e+00),
          30.041069500000, 0.00, "R2"},
         {/*275.555*/ DoodsonConstituent(2, 2, 0, 0, 0, 0, +0.0e+00),
          30.082140100000, 0.00, "K2"},
          /* extraced from [3] Table 5.1 */
         {/*245.655*/ DoodsonConstituent(2, -1, 0, 1, 0, 0, +0.0e+00), 
          28.439729500000, 0.00, "N2"},
          /* extraced from [3] Table 5.1 */
         {/*265.455*/ DoodsonConstituent(2, 1, 0, -1, 0, 0, +2.0e+00), 
          28.439729500000, 0.00, "L2"},
         {/*381.555*/ DoodsonConstituent(3, 3, -4, 0, 0, 0, +0.0e+00),
          44.958935590000, 0.00, "T3"},
         {/*382.555*/ DoodsonConstituent(3, 3, -3, 0, 0, 0, +2.0e+00),
          45.000004230000, 0.00, "S3"},
         {/*383.555*/ DoodsonConstituent(3, 3, -2, 0, 0, 0, +0.0e+00),
          45.041072870000, 0.00, "R3"},
         {/*491.555*/ DoodsonConstituent(4, 4, -4, 0, 0, 0, +0.0e+00),
          60.000005640000, 0.00, "S4"},
         {/*5100.555*/ DoodsonConstituent(5, 5, -5, 0, 0, 0, +0.0e+00),
          75.000007050000, 0.00, "S5"},
         {/*611-1.555*/ DoodsonConstituent(6, 6, -6, 0, 0, 0, +0.0e+00),
          90.000008460000, 0.00, "S6"}}}; /* AtmosphericTidalHarmonics */

constexpr static std::array<TidalConstituentArrayEntry, 34>
  OceanicTidalHarmonics = {
{{/*055.565*/DoodsonConstituent(0,0,0,0,1,0,2e0),  0e0, 0e0, "Omega1"},/* [4] and [5], EOT11 */
 {/*055.575*/DoodsonConstituent(0,0,0,0,2,0,0),    0e0, 0e0, "Omega2"},/* [4] and [5], EOT11 */
 {/*056.554*/DoodsonConstituent(0,0,1,0,0,-1,0),   0e0, 0e0, "Sa"    },/* [4] and [5], EOT11 */
 {/*057.555*/DoodsonConstituent(0,0,2,0,0,0,0),    0e0, 0e0, "Ssa"   },/* [4] and [5], EOT11 */
 {/*065.455*/DoodsonConstituent(0,1,0,-1,0,0,0),   0e0, 0e0, "Mm"    },/* [4] and [5], EOT11 */
 {/*075.555*/DoodsonConstituent(0,2,0,0,0,0,0),    0e0, 0e0, "Mf"    },/* [4] and [5], EOT11 */
 {/*085.455*/DoodsonConstituent(0,3,0,-1,0,0,0),   0e0, 0e0, "Mtm"   },/* [4] and [5], EOT11 */
 {/*093.555*/DoodsonConstituent(0,4,-2,0,0,0,0),   0e0, 0e0, "Msqm"  },/* [4] and [5], EOT11 */
 
 {/*135.655*/DoodsonConstituent(1,-2,0,1,0,0,-1e0), 0e0, 0e0, "Q1"   },/* [4] and [5], EOT11 */
 {/*145.555*/DoodsonConstituent(1,-1,0,0,0,0,-1e0), 0e0, 0e0, "O1"   },/* [4] and [5], EOT11 */
 {/*163.555*/DoodsonConstituent(1,1,-2,0,0,0,-1e0), 0e0, 0e0, "P1"   },/* [4] and [5], EOT11 */
 {/*165.555*/DoodsonConstituent(1,1,0,0,0,0,+1e0),  0e0, 0e0, "K1"   },/* [4] and [5], EOT11 */
 {/*164.555*/DoodsonConstituent(1,1,-1,0,0,0,+2e0), 0e0, +0e+00, "S1"},/* [5] */
 {/*175.455*/DoodsonConstituent(1,2,0,-1,0,0,+1e0), 0e0, +0e+00, "J1"},/* [5] */

 
 {/*235.755*/DoodsonConstituent(2,-2,0,2,0,0,0),   0e0, 0e0, "2N2"   },/* [4] and [5], EOT11 */
 {/*245.655*/DoodsonConstituent(2,-1,0,1,0,0,0),   0e0, 0e0, "N2"    },/* [4] and [5], EOT11 */
 {/*255.555*/DoodsonConstituent(2,0,0,0,0,0,0),    0e0, 0e0, "M2"    },/* [4] and [5], EOT11 */
 {/*273.555*/DoodsonConstituent(2,2,-2,0,0,0,0),   0e0, 0e0, "S2"    },/* [4] and [5], EOT11 */
 {/*275.555*/DoodsonConstituent(2,2,0,0,0,0,0),    0e0, 0e0, "K2"    },/* [4] and [5], EOT11 */
 {/*272.556*/DoodsonConstituent(2,2,-3,0,0,1,0),   0e0, +0e+00,   "T2" },/* [5] */
 {/*227.655*/DoodsonConstituent(2,-3,2,1,0,0,0),   0e0, +0e+00,   "Epsilon2"},/* [5] */
 {/*237.555*/DoodsonConstituent(2,-2,2,0,0,0,0),   0e0, +0e+00,   "Mu2"}, /* [5] */
 {/*247.455*/DoodsonConstituent(2,-1,2,-1,0,0,0),  0e0, +0e+00,   "Nu2"},/* [5] */
 {/*263.655*/DoodsonConstituent(2,1,-2,1,0,0,+2e0), 0e0, +0e+00,  "La2"},/* [5] */
 {/*265.455*/DoodsonConstituent(2,1,0,-1,0,0,+2e0), 0e0, +0e+00,  "L2"},/* [5] */
 {/*274.554*/DoodsonConstituent(2,2,-1,0,0,-1,+2e0), 0e0, +0e+00, "R2"},/* [5] */

 {/*355.555*/DoodsonConstituent(3,0,0,0,0,0,+2e0), 0e0, +0e+00, "M3"},/* [5] */
 
 {/*455.555*/DoodsonConstituent(4,0,0,0,0,0,0),  0e0, 0e0,    "M4"    },/* [4] and [5], EOT11 */
 {/*435.755*/DoodsonConstituent(4,-2,0,2,0,0,0), 0e0, +0e+00, "N4"},/* [5] */
 {/*445.655*/DoodsonConstituent(4,-1,0,1,0,0,0), 0e0, +0e+00, "Mn4"},/* [5] */
 {/*473.555*/DoodsonConstituent(4,2,-2,0,0,0,0), 0e0, +0e+00, "Ms4"},/* [5] */
 {/*491.555*/DoodsonConstituent(4,4,-4,0,0,0,0), 0e0, +0e+00, "S4"},/* [5] */

 {/*655.555*/DoodsonConstituent(6,0,0,0,0,0,0),  0e0, +0e+00, "M6"},/* [5] */
 
 {/*855.555*/DoodsonConstituent(8,0,0,0,0,0,0),  0e0, +0e+00, "M8"},/* [5] */
}};

const TidalConstituentArrayEntry *
match_ocean_tide_wave(const DoodsonConstituent &wave) noexcept;
const TidalConstituentArrayEntry *
match_ocean_tide_wave(const char *wave) noexcept;
} /* namespace detail */

class TidalWave {
  static constexpr const int max_chars = 15; /* +1 for '\0' */

  DoodsonConstituent _d;
  double _per;            /* period in [deg/hour] */
  double _hf;             /* height */
  char _n[max_chars + 1]; /* name */

public:
  const DoodsonConstituent &doodson() const noexcept { return _d; }
  DoodsonConstituent &doodson() noexcept { return _d; }
  double period() const noexcept { return _per; }
  double &period() noexcept { return _per; }
  double height() const noexcept { return _hf; }
  double &height() noexcept { return _hf; }
  const char *name() const noexcept { return _n; }

  TidalWave() noexcept {};

  TidalWave(const detail::TidalConstituentArrayEntry &entry) noexcept
      : _d(entry._d), _per(entry._per), _hf(entry._hf) {
    std::strcpy(_n, entry._n);
  }

  TidalWave(const DoodsonConstituent &d, double per = 0e0, double hf = 0e0,
            const char *name = nullptr) noexcept
      : _d(d), _per(per), _hf(hf) {
    std::memset(_n, '\0', max_chars + 1);
    if (name) {
      const int sz = std::strlen(name);
      std::memcpy(_n, name, sz > max_chars ? max_chars : sz);
    }
  }
}; /* TidalWave */

/** Given a tidal constituent name, return its details */
const detail::TidalConstituentArrayEntry *
find_wave_entry(const char *name) noexcept;
TidalWave get_wave(const char *name);
} /* namespace dso */
#endif
