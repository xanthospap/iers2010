/** @file
 * Atmospheric Tide Model & computation for Geodesy.
 */

#ifndef __DSO_SOLID_EARTH_TIDE_HPP__
#define __DSO_SOLID_EARTH_TIDE_HPP__

#include "doodson.hpp"
#include "aod1b.hpp"
#include "stokes_coefficients.hpp"
#include "datetime/calendar.hpp"
#include "geodesy/geodesy.hpp"
#include "eigen3/Eigen/Eigen"
#include <vector>

namespace dso {

namespace detail {
class AtmosphericTidalWave {
  TidalConstituentsArrayEntry mdentry;
  StokesCoeffs mCosCs;
  StokesCoeffs mSinCs;
}; /* AtmosphericTidalWave */
} /* namespace detail */

class AtmosphericTides {
private:
  std::vector<detail::AtmosphericTidalWave> mwaves;
public:
  int append_wave(const char *aod1b_fn) noexcept {
    AtmosphericTidalWave newWave;

    return Aod1bIn::get_tidal_wave_coeffs(aod1b_fn, ,
                                   newWave.mCosCs, newWave.msinCs, Aod1bIn *aod1b,
                                   int max_degree = -1,
                                   int max_order = -1) noexcept;
  }
}; /* AtmosphericTides */

}/* namespace dso */

#endif
