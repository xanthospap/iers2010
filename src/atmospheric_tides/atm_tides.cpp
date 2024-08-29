#include "atmospheric_tides.hpp"
#include <cstdio>
#include <stdexcept>

std::vector<dso::detail::AtmosphericTidalWave>::iterator
dso::AtmosphericTide::append_wave(
    const dso::TidalWave &wave, int max_degree,
    int max_order) {
  if (auto it = find_tidal_wave(wave.doodson()); it == mwaves.end()) {
    /* add new constituent */
    mwaves.emplace_back(
        dso::detail::AtmosphericTidalWave(wave, max_degree, max_order));
  } else {
    /* this doodson number (wave) already exists */
    char buf[8] = "\0";
    fprintf(stderr,
            "[ERROR] Tidal wave with doodson number %s already exists; failed "
            "to append to AtmosphericTide instance (traceback: %s)\n",
            buf, __func__);
    throw std::runtime_error("[ERROR] Failed pushing back tidal wave\n");
  }
  return mwaves.end() - 1;
}
