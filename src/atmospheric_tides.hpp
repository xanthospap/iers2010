/** @file
 * Atmospheric Tide Model & computation for Geodesy.
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

namespace detail {
class AtmosphericTidalWave {
private:
  TidalWave mwave;
  StokesCoeffs mCosCs;
  StokesCoeffs mSinCs;

public:
  const StokesCoeffs &stokes_sin() const noexcept { return mSinCs; }
  StokesCoeffs &stokes_sin() noexcept { return mSinCs; }
  const StokesCoeffs &stokes_cos() const noexcept { return mCosCs; }
  StokesCoeffs &stokes_cos() noexcept { return mCosCs; }
  TidalWave wave() const noexcept { return mwave; }
  TidalWave &wave() noexcept { return mwave; }
  
  AtmosphericTidalWave(const TidalConstituentsArrayEntry *wave, double Gm,
                       double Re, int max_degree, int max_order) noexcept
      : mwave(*wave), mCosCs(max_degree, max_order, Gm, Re),
        mSinCs(max_degree, max_order, Gm, Re){};

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
