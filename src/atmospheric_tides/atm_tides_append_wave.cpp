#include "atmospheric_tides.hpp"
#include <stdexcept>

int dso::AtmosphericTides::append_wave(const char *aod1b_fn, int max_degree,
                                       int max_order) noexcept {
  /* max degree and order collected */
  int mdeg_collected = 0;
  int mord_collected = 0;

  try {
    /* constructor may throw */
    dso::Aod1bIn aod1b(aod1b_fn);

    /* set degree and order of coeffs */
    if (max_degree < 0)
      max_degree = aod1b.max_degree();
    if (max_order < 0)
      max_order = aod1b.max_degree();

    /* tidal wave in AOD1B file not resolved! */
    if (!aod1b.tidal_wave())
      return 1;

    /* construct new instance */
    dso::detail::AtmosphericTidalWave newWave(
        aod1b.tidal_wave(), aod1b.GM(), aod1b.Re(), max_degree, max_order);

    /* parse coefficients */
    if (aod1b.get_tidal_wave_coeffs(newWave.mCosCs, newWave.mSinCs, max_degree,
                                    max_order)) {
      return 1;
    }

    /* push back the new wave */
    mwaves.emplace_back(newWave);

    /* set degree and order collected */
    mdeg_collected = newWave.mCosCs.max_degree();
    mord_collected = newWave.mCosCs.max_order();

  } catch (std::runtime_error &) {
    return 1;
  }

  /* do we need to re-size the instance's stokes coeffs size? */
  if (mdeg_collected > mcs.max_degree() || mord_collected > mcs.max_order())
    mcs.cresize(mdeg_collected, mord_collected);
  if (mdeg_collected > mcs.max_degree() || mord_collected > mcs.max_order())
    mcs.cresize(mdeg_collected, mord_collected);

  return 0;
}
