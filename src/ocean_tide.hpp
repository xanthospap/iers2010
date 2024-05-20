/** @file
 */

#ifndef __DSO_OCEAN_TIDE_HPP__
#define __DSO_OCEAN_TIDE_HPP__

#include "aod1b.hpp"
#include "datetime/calendar.hpp"
#include "doodson.hpp"
#include "geodesy/geodesy.hpp"
#include "stokes_coefficients.hpp"
#include <algorithm>
#include <vector>

namespace dso {

namespace detail {
class OceanicTidalWave {
public:
private:
  DoodsonConstituent mdoodson;
  StokesCoeffs mCosCs;
  StokesCoeffs mSinCs;

public:
  OceanicTidalWave(const DoodsonConstituent &wave, double Gm, double Re,
                   int max_degree, int max_order) noexcept
      : mdoodson(wave), mCosCs(max_degree, max_order, Gm, Re),
        mSinCs(max_degree, max_order, Gm, Re){};

  StokesCoeffs &cos_coeffs() noexcept { return mCosCs; }
  StokesCoeffs &sin_coeffs() noexcept { return mSinCs; }
  DoodsonConstituent &doodson() noexcept { return mdoodson; }
  const StokesCoeffs &cos_coeffs() const noexcept { return mCosCs; }
  const StokesCoeffs &sin_coeffs() const noexcept { return mSinCs; }
  const DoodsonConstituent &doodson() const noexcept { return mdoodson; }
}; /* OceanicTidalWave */
} /* namespace detail */

class OceanTide {
private:
  std::vector<detail::OceanicTidalWave> mwaves;
  StokesCoeffs mcs;

public:
  auto find_wave(const DoodsonConstituent &d) noexcept {
    return std::find_if(
        mwaves.begin(), mwaves.end(),
        [&](const detail::OceanicTidalWave &w) { return w.doodson() == d; });
  }

  const std::vector<detail::OceanicTidalWave> &waves() const noexcept {
    return mwaves;
  }
  
  std::vector<detail::OceanicTidalWave> &waves() noexcept {
    return mwaves;
  }

  void reserve(int num_waves) noexcept {mwaves.reserve(num_waves);}

  void append_wave(detail::OceanicTidalWave &&wave) noexcept {
    mwaves.emplace_back(wave);
  }

  void append_wave(const DoodsonConstituent &doodson, int max_degree,
                   int max_order) noexcept {
    mwaves.emplace_back(detail::OceanicTidalWave(
        doodson, iers2010::GMe, iers2010::Re, max_degree, max_order));
  }

  int stokes_coeffs(int max_degree, int max_order, const MjdEpoch &mjdtt,
                    const MjdEpoch &mjdut1,
                    const double *const delaunay_args) noexcept;

}; /* OceanTide */

OceanTide initFes2014bFromIcgem(const char *dir, const char *fn_generic_name,
                                int max_degree = 180, int max_order = 180);

} /* namespace dso */

#endif
