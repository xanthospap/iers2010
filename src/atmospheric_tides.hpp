/** @file
 * Atmospheric Tide Model & computation for Geodesy.
 *
 * References:
 *
 * [1] Sulzbach, Roman; Balidakis, Kyriakos; Dobslaw, Henryk; Thomas,
 * Maik (2022): TiME22: Periodic Disturbances of the Terrestrial Gravity
 * Potential Induced by Oceanic and Atmospheric Tides. GFZ Data Services.
 * https://doi.org/10.5880/GFZ.1.3.2022.006
 *
 * [2] Henryk Dobslaw, Inga Bergmann-Wolf, Robert Dill, Lea Poropat,
 * Frank Flechtner (2017): GRACE 327-750, Gravity Recovery and Climate
 * Experiment, Product Description Document for AOD1B Release 06 (Rev. 6.1,
 * October 19, 2017), GFZ German Research Centre for Geosciences Department 1:
 * Geodesy
 *
 * [3] Balidakis, K., Sulzbach, R., Shihora, L., Dahle, C., Dill, R., &
 * Dobslaw, H. (2022). Atmospheric contributions to global ocean tides for
 * satellite gravimetry. Journal of Advances in Modeling Earth Systems, 14,
 * e2022MS003193. https://doi.org/10.1029/2022MS003193
 *
 * [4] Lasser, M., Meyer, U., Jäggi, A., Mayer-Gürr, T., Kvas, A., Neumayer,
 * K. H., Dahle, C., Flechtner, F., Lemoine, J.-M., Koch, I., Weigelt, M.,
 * and Flury, J.: Benchmark data for verifying background model
 * implementations in orbit and gravity field determination software,
 * Adv. Geosci., 55, 1–11, https://doi.org/10.5194/adgeo-55-1-2020, 2020.
 *
 * [5] GFZ-Potsdam Information System And Data Center, Global Earth Science
 * Data, ftp://isdcftp.gfz-potsdam.de/grace/Level-1B/GFZ/AOD
 */

#ifndef __DSO_ATMOSPHERIC_TIDE_HPP__
#define __DSO_ATMOSPHERIC_TIDE_HPP__

#include "aod1b.hpp"
#include "tide_atlas.hpp"
#include "datetime/calendar.hpp"
#include <vector>

namespace dso {

/** @class AtmosphericTide
 *
 * A class to hold an Atmospheric Tide Atlas, i.e. a model to assist computing
 * the gravitational effect of atmospheric tides.
 * A model usually includes a series of tidal constituents (waves), each of
 * which has an individual contribution to the final effect.
 * At any given time, the gravitational effect can be represented as
 * corrections to the Stokes coefficients of a spherical harmonics expansion.
 * Any tidal constituent included in the model, should be an instance of the
 * AtmosphericTidalWave, so that it holds in-phase and quadrature Stokes
 * coeffcients.
 * For the computation of the (corrections to the) Stokes coefficients, we
 * follow Eq. (28) of [4].
 */
class AtmosphericTide {
private:
  /* The tide atlas */
  TideAtlas matlas;
  /* Stokes coeffs to hold (if needed) the accumulated effect of all waves */
  StokesCoeffs mcs;

public:
  const TideAtlas &atlas() const noexcept {return matlas;}
  TideAtlas &atlas() noexcept {return matlas;}

  AtmosphericTide(const TideAtlas &atlas, const char *name = nullptr)
      : matlas(atlas), mcs(atlas.max_atlas_degree()) {
    int max_order;
    int max_degree = matlas.max_atlas_degree(max_order);
    if (max_order < max_degree)
      mcs.resize(max_degree, max_order);
    if (name)
      std::strcpy(matlas.name(), name);
  }

  /** @brief Append a tidal wave using the AOD1B format.
   *
   * This function is normally used to append a specific wave, via an AOD1B
   * RLXX product file. It will read and parse the file, resolve the
   * constituent (AOD1B files only contain the name of the wave, not the
   * Doodson number) and store the in-phase and quadrature (Stokes) coeffs.
   *
   * Note that AOD1B only contain the tidal constituent name (e.g. "K1") not
   * the Doodson number. To get the Doodson number (which we need), the name
   * recorded in the product file should match to an entry in the
   * detail::AtmosphericTidalHarmonics table (see file doodson.hpp). If the
   * name is not matched, an error will occur (non-zero return code).
   *
   * @param[in] aod1b_fn An AOD1B product file of some atmospheric tidal
   *            constituent. Such files are available e.g. at [5]
   * @param[in] max_degree Max degree of in-phase and quadrature (Stokes)
   *            coeffs to read/parse/store for the new wave.
   * @param[in] max_order  Max order of in-phase and quadrature (Stokes)
   *            coeffs to read/parse/store for the new wave.
   * @return Anything other than zero denotes an error.
   */
  int append_wave(const char *aod1b_fn, int max_degree, int max_order) noexcept;

  /** @brief Compute accumulated Stokes coefficients.
   *
   * Compute the accumulated Stokes coefficients using all available waves
   * included in the instance. The computation follows Eq. (28) of [4].
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

  /** @brief Get the instance's Stokes coefficients */
  const StokesCoeffs &stokes_coeffs() const noexcept { return mcs; }

  /** @brief Get the instance's Stokes coefficients */
  StokesCoeffs &stokes_coeffs() noexcept { return mcs; }
}; /* AtmosphericTides */

} /* namespace dso */

#endif
