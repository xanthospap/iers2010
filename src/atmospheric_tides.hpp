/** @file
 * Atmospheric Tide Model & computation for Geodesy.
 */

#ifndef __DSO_SOLID_EARTH_TIDE_HPP__
#define __DSO_SOLID_EARTH_TIDE_HPP__

#include "aod1b.hpp"
#include "datetime/calendar.hpp"
#include "doodson.hpp"
#include "geodesy/geodesy.hpp"
#include "stokes_coefficients.hpp"
#include <vector>

namespace dso {

namespace detail {
class AtmosphericTidalWave {
public:
  TidalConstituentsArrayEntry mdentry;
  StokesCoeffs mCosCs;
  StokesCoeffs mSinCs;
  AtmosphericTidalWave(const TidalConstituentsArrayEntry *wave, double Gm,
                       double Re, int max_degree, int max_order) noexcept
      : mdentry(*wave), mCosCs(max_degree, max_order, Gm, Re),
        mSinCs(max_degree, max_order, Gm, Re){};
}; /* AtmosphericTidalWave */
} /* namespace detail */

class AtmosphericTides {
private:
  std::vector<detail::AtmosphericTidalWave> mwaves;
  StokesCoeffs mcs;

public:
  int append_wave(const char *aod1b_fn, int max_degree, int max_order) noexcept;

  int stokes_coeffs(const MjdEpoch &mjdtt, const MjdEpoch &mjdut1,
                    const double *const delaunay_args) noexcept;

  const StokesCoeffs &stokes_coeffs() const noexcept { return mcs; }
  StokesCoeffs &stokes_coeffs() noexcept { return mcs; }
}; /* AtmosphericTides */

} /* namespace dso */

#endif
