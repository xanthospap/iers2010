#include "atmospheric_tides.hpp"
#include <stdexcept>

int dso::AtmosphericTide::append_wave(const char *aod1b_fn, int max_degree,
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

    /* construct new instance (may throw if not resolved) */
    dso::TidalWave newWave(aod1b.tidal_wave());

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
  if (mdeg_collected > mcs.max_degree() || mord_collected > mcs.max_order()) {
    mcs.cresize(mdeg_collected, mord_collected);
  }
  if (mdeg_collected > mcs.max_degree() || mord_collected > mcs.max_order()) {
    mcs.cresize(mdeg_collected, mord_collected);
  }

  //int sz = mwaves.size();
  //printf("Collected wave: %s with factor: %.1f, d/o: %d/%d\n",
  //       mwaves[sz - 1].mdentry._n, mwaves[sz - 1].mdentry._d.pifactor(),
  //       mwaves[sz - 1].mCosCs.max_degree(),
  //       mwaves[sz - 1].mCosCs.max_order());

  return 0;
}
