/** @file
 */

#ifndef __DSO_TIDE_ATLAS_HPP__
#define __DSO_TIDE_ATLAS_HPP__

#include "doodson.hpp"
#include "groops_tide_model.hpp"
#include "stokes_coefficients.hpp"
#include <vector>

namespace dso {

/** @class TidalConstituent aka a tidal line.
 *
 * Why do we need this class?
 * Well, each ocean and atmospheric model, can be given as an atlas, i.e. a
 * list of major tidal lines/constituents, complemented by sin- and cos-
 * stokes coefficients for a shperical harmonics expansion. Adding the
 * effects of the tidal lines, will give the total effect of the model
 * (possibly including admittance).
 *
 * Each TidalConstituent represents a major tidal line of a tide atlas. It
 * should hold information on the wave itself (i.e. doodson number, period,
 * etc) encoded into a TidalWave instance, and ofcourse it should also hold
 * the in- and out-of- phase stokes coefficients.
 */
class TidalConstituent {
private:
  /* tidal wave information */
  TidalWave mwave;
  /* in-phase Stokes coefficients, i.e. Cnm_cos and Snm_cos */
  StokesCoeffs mCosCs;
  /* quadrature Stokes coefficients, i.e. Cnm_sin and Snm_sin */
  StokesCoeffs mSinCs;

public:
  const StokesCoeffs &stokes_sin() const noexcept { return mSinCs; }
  StokesCoeffs &stokes_sin() noexcept { return mSinCs; }
  const StokesCoeffs &stokes_cos() const noexcept { return mCosCs; }
  StokesCoeffs &stokes_cos() noexcept { return mCosCs; }
  TidalWave wave() const noexcept { return mwave; }
  TidalWave &wave() noexcept { return mwave; }
  const DoodsonConstituent &doodson() const { return mwave.doodson(); }
  const char *name() const { return mwave.name(); }

  /* Constructor given a TidalConstituentsArrayEntry and the maximum degree
   * and order of the relevant Stokes coefficients.
   */
  TidalConstituent(const detail::TidalConstituentArrayEntry *wave,
                   int max_degree, int max_order) noexcept;

  /* Constructor given a TidalWave and the maximum degree and order of the
   * relevant Stokes coefficients.
   */
  TidalConstituent(const TidalWave &wave, int max_degree,
                   int max_order) noexcept;

}; /* TidalConstituent */

/** @class TidalAtlas
 */
class TideAtlas {
public:
  static constexpr const int NAME_MAX_CHARS = 64;

private:
  /* list of tidal constituents (waves) included in the model */
  std::vector<TidalConstituent> mwaves;
  /* name of the model */
  char mname[NAME_MAX_CHARS] = {'\0'};

public:
  int max_atlas_degree(int &max_order) const noexcept {
    int max_degree = 0;
    max_order = 0;
    for (const auto &w : mwaves) {
      max_degree = std::max(max_degree, w.stokes_cos().max_degree());
      max_order = std::max(max_order, w.stokes_cos().max_order());
    }
    return max_degree;
  }
  int max_atlas_degree() const noexcept {
    int max_degree = 0;
    for (const auto &w : mwaves) {
      max_degree = std::max(max_degree, w.stokes_cos().max_degree());
    }
    return max_degree;
  }

  /* @brief Search the list of tidal constituent for a specific wave.
   *
   * Given a Doodson number, search through the list of constituents included
   * in the instance, and return the one matching the Doodson number.
   * For more details on "matching", see the DoodsonConstituent class
   * operators.
   *
   * @param[in] doodson The DoodsonConstituent to search for
   * @return If the given DoodsonConstituent is matched, return an iterator
   *         to the instance's mwaves vector, that points to the matching wave.
   *         If not matched, then return mwaves.end().
   */
  std::vector<TidalConstituent>::const_iterator
  find_tidal_wave(const DoodsonConstituent &doodson) const noexcept;
  std::vector<TidalConstituent>::iterator
  find_tidal_wave(const DoodsonConstituent &doodson) noexcept;

  /** Return the vector of waves, i.e. mwaves */
  const std::vector<TidalConstituent> &waves() const noexcept { return mwaves; }

  /** Return the vector of waves, i.e. mwaves */
  std::vector<TidalConstituent> &waves() noexcept { return mwaves; }

  const char *name() const { return mname; }
  char *name() { return mname; }

  /** @brief Append a tidal wave.
   *
   * Append a tidal wave using an TidalWave instance to describe it. The new
   * tidal wave will have the specified degree and order of (in-phase and
   * quadrature) Stokes coefficients, which should be set later on.
   *
   * If the wave to be added has a DoodonsonConstituent that is already
   * included in the mwaves, it is considered as duplicate and an exception
   * will be thrown.
   *
   * @param[in] wave A TidalWave instance that describes the wave to be added
   * @param[in] max_degree Max degree of the in-phase and quadrature Stokes
   *                 coefficients
   * @param[in] max_order Max order of the in-phase and quadrature Stokes
   *                 coefficients
   * @return On success (i.e. the wave was added) the function will return an
   * iterator to the newly added wave within the mwaves vector. If the wave is
   * already included in the mwaves vector, an exception will be thrown.
   */
  std::vector<TidalConstituent>::iterator
  append_wave(const TidalConstituent &wave);

  std::vector<TidalConstituent>::iterator
  append_wave(const TidalWave &wave, int max_degree, int max_order);

}; /* TideAtlas */

TideAtlas groops_atlas(const char *fn, const char *dir, int max_degree = -1,
                       int max_order = -1);
TideAtlas groops_atlas(const char *file_001, const char *file_002,
                       const char *file_003, const char *dir,
                       GroopsTideModel &mdl, int max_degree = -1,
                       int max_order = -1);

} /* namespace dso */

#endif
