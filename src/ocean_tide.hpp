/** @file
 *  Ocean Tide Model & computation for Geodesy.
 *
 *  References:
 *
 * [1] Lasser, M., Meyer, U., Jäggi, A., Mayer-Gürr, T., Kvas, A., Neumayer,
 * K. H., Dahle, C., Flechtner, F., Lemoine, J.-M., Koch, I., Weigelt, M.,
 * and Flury, J.: Benchmark data for verifying background model
 * implementations in orbit and gravity field determination software,
 * Adv. Geosci., 55, 1–11, https://doi.org/10.5194/adgeo-55-1-2020, 2020.
 */

#ifndef __DSO_OCEAN_TIDE_HPP__
#define __DSO_OCEAN_TIDE_HPP__

#include "aod1b.hpp"
#include "datetime/calendar.hpp"
#include "doodson.hpp"
#include "geodesy/transformations.hpp"
#include "stokes_coefficients.hpp"
#include <algorithm>
#include <vector>

namespace dso {

namespace detail {

/** @class OceanicTidalWave
 *
 * A tidal constituent (wave) to represent contribution from an individual
 * wave. A wave has a TidalWave component that holds its characteristics
 * (e.g. Doodson number) as well as Stokes coefficients for in-phase (Cnm_cos
 * and Snm_cos) and quadrature (Cnm_sin and Snm_sin) components.
 *
 * For example, an instance of OceanicTidalWave could be used to compute
 * the potential (via spherical harmonics) due to oceanic tidal loading
 * owing to the M2 tidal wave.
 *
 * Normally, an Oceanic Loading model includes multiple waves, i.e.
 * multiple OceanicTidalWave's constituting a tidal atlas.
 *
 * References: [1], [2]
 */
class OceanicTidalWave {
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

  /* Constructor given a TidalConstituentsArrayEntry and the maximum degree
   * and order of the relevant Stokes coefficients.
   */
  OceanicTidalWave(const TidalConstituentArrayEntry *wave, int max_degree,
                       int max_order) noexcept
      : mwave(*wave), mCosCs(max_degree, max_order),
        mSinCs(max_degree, max_order) {};

  /* Constructor given a TidalWave and the maximum degree and order of the
   * relevant Stokes coefficients.
   */
  OceanicTidalWave(const TidalWave &wave, int max_degree,
                       int max_order) noexcept
      : mwave(wave), mCosCs(max_degree, max_order),
        mSinCs(max_degree, max_order) {};
  OceanicTidalWave(const DoodsonConstituent &wave, double Gm, double Re,
                   int max_degree, int max_order) noexcept
      : mwave(wave), mCosCs(max_degree, max_order, Gm, Re),
        mSinCs(max_degree, max_order, Gm, Re){};

}; /* OceanicTidalWave */
} /* namespace detail */

/** @class OceanTide
 *
 * A class to hold an Ocean Tide Atlas, i.e. a model to assist computing
 * the gravitational effect of oceanic tides.
 * A model usually includes a series of tidal constituents (waves), each of
 * which has an individual contribution to the final effect.
 * At any given time, the gravitational effect can be represented as
 * corrections to the Stokes coefficients of a spherical harmonics expansion.
 * Any tidal constituent included in the model, should be an instance of the
 * OceanicTidalWave, so that it holds in-phase and quadrature Stokes
 * coeffcients.
 * For the computation of the (corrections to the) Stokes coefficients, we
 * follow Eq. (17) of [1]
 */
class OceanTide {
  static constexpr const int NAME_MAX_CHARS = 64;
private:
  /* vector of individual waves */
  std::vector<detail::OceanicTidalWave> mwaves;
  /* Stokes coeffs to hold (if needed) the accumulated effect of all waves */
  StokesCoeffs mcs;
  /* model name */
  char mname[NAME_MAX_CHARS] = {'\0'};

public:
  
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
  auto find_tidal_wave(const DoodsonConstituent &d) noexcept {
    return std::find_if(mwaves.begin(), mwaves.end(),
                        [&](const detail::OceanicTidalWave &w) {
                          return w.wave().doodson() == d;
                        });
  }

  const std::vector<detail::OceanicTidalWave> &waves() const noexcept {
    return mwaves;
  }
  
  /** Return the vector of waves, i.e. mwaves */
  std::vector<detail::OceanicTidalWave> &waves() noexcept {
    return mwaves;
  }

  /* reserve number of individual tidal waves in atlas */
  void reserve(int num_waves) noexcept {mwaves.reserve(num_waves);}

  /** @brief Append a tidal wave.
   *
   * Append a tidal wave using an OceanicTidalWave instance to describe it. 
   *
   * If the wave to be added has a DoodonsonConstituent that is already
   * included in the mwaves, it is considered as duplicate and an exception
   * will be thrown.
   *
   * @param[in] wave An OceanicTidalWave instance that describes the wave to 
   *                 be added
   * @return On success (i.e. the wave was added) the function will return an
   * iterator to the newly added wave within the mwaves vector. If the wave is
   * already included in the mwaves vector, an exception will be thrown.
   */
  std::vector<detail::OceanicTidalWave>::iterator
  append_wave(detail::OceanicTidalWave &&wave);

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
  std::vector<detail::OceanicTidalWave>::iterator
  append_wave(const TidalWave &wave, int max_degree, int max_order);

  /** @brief Compute accumulated Stokes coefficients.
   *
   * Compute the accumulated Stokes coefficients using all available waves
   * included in the instance. The computation follows Eq. (17) of [1].
   * The computed Stokes coefficients are available in hte instance's mcs
   * member variable.
   *
   * @param[in] mjdtt Epoch for computation in TT
   * @param[in] mjdtt Epoch for computation in UT1
   * @param[in] delaunay_args An array containing Delaunay (i.e. fundamental)
   *            arguments, in the order:
   *            l  : Mean anomaly of the Moon [rad]
   *            lp : Mean anomaly of the Sun [rad]
   *            f  : L - OM [rad]
   *            d  : Mean elongation of the Moon from the Sun [rad]
   *            om : Mean longitude of the ascending node of the Moon [rad]
   * @return Always zero. Note that the computed Stokes coeffs are stored in
   *            the instance's mcs member variable.
   */
  int stokes_coeffs(const MjdEpoch &mjdtt, const MjdEpoch &mjdut1,
                    const double *const delaunay_args) noexcept;
  
  const StokesCoeffs &stokes_coeffs() const noexcept {return mcs;}
  StokesCoeffs &stokes_coeffs() noexcept {return mcs;}
  
  void resize_stokes_ceoffs(int max_degree, int max_order) noexcept {
    mcs.resize(max_degree, max_order);
  }

  char *name() noexcept {return mname;}
  const char *name() const noexcept {return mname;}

}; /* OceanTide */

OceanTide initFes2014bFromIcgem(const char *dir, const char *fn_generic_name,
                                int max_degree = -1, int max_order = -1);

} /* namespace dso */

#endif
