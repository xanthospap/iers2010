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
 */

#ifndef __DSO_SOLID_EARTH_TIDE_HPP__
#define __DSO_SOLID_EARTH_TIDE_HPP__

#include "aod1b.hpp"
#include "datetime/calendar.hpp"
#include "doodson.hpp"
#include "geodesy/transformations.hpp"
#include "stokes_coefficients.hpp"
#include <vector>

namespace dso {

/** @class AtmosphericTidalWave
 *
 * A tidal constituent (wave) to represent contribution from an individual 
 * wave. A wave has a TidalWave component that holds its characteristics 
 * (e.g. Doodson number) as well as Stokes coefficients for in-phase (Cnm_cos 
 * and Snm_cos) and quadrature (Cnm_sin and Snm_sin) components.
 *
 * References: [1], [2]
 */
namespace detail {
class AtmosphericTidalWave {
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
  AtmosphericTidalWave(const TidalConstituentArrayEntry *wave, int max_degree,
                       int max_order) noexcept
      : mwave(*wave), mCosCs(max_degree, max_order),
        mSinCs(max_degree, max_order) {};

  /* Constructor given a TidalWave and the maximum degree and order of the 
   * relevant Stokes coefficients.
   */
  AtmosphericTidalWave(const TidalWave &wave, int max_degree,
                       int max_order) noexcept
      : mwave(wave), mCosCs(max_degree, max_order),
        mSinCs(max_degree, max_order) {};

}; /* AtmosphericTidalWave */
} /* namespace detail */

class AtmosphericTide {
private:
  std::vector<detail::AtmosphericTidalWave> mwaves;
  StokesCoeffs mcs;

public:
  auto find_tidal_wave(const DoodsonConstituent &doodson) const noexcept {
    return std::find_if(mwaves.begin(), mwaves.end(),
                        [=](const detail::AtmosphericTidalWave &w) {
                          return w.wave().doodson() == doodson;
                        });
  }

  const auto wave_vector() const noexcept { return mwaves; }
  auto wave_vector() noexcept { return mwaves; }

  int append_wave(const char *aod1b_fn, int max_degree, int max_order) noexcept;

  std::vector<detail::AtmosphericTidalWave>::iterator
  append_wave(const TidalWave &wave, int max_degree,
              int max_order);

  int stokes_coeffs(const MjdEpoch &mjdtt, const MjdEpoch &mjdut1,
                    const double *const delaunay_args) noexcept;

  const StokesCoeffs &stokes_coeffs() const noexcept { return mcs; }
  StokesCoeffs &stokes_coeffs() noexcept { return mcs; }
}; /* AtmosphericTides */

} /* namespace dso */

#endif
