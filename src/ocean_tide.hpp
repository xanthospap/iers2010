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
#include "tide_atlas.hpp"

namespace dso {

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
private:
  /* the tidal atlas */
  TideAtlas matlas;
  /* Stokes coeffs to hold (if needed) the accumulated effect of all waves */
  StokesCoeffs mcs;

public:
  const TideAtlas &atlas() const noexcept {return matlas;}
  TideAtlas &atlas() noexcept {return matlas;}
  
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

}; /* OceanTide */

OceanTide initFes2014bFromIcgem(const char *dir, const char *fn_generic_name,
                                int max_degree = -1, int max_order = -1);

} /* namespace dso */

#endif
