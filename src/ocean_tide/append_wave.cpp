#include "ocean_tide.hpp"

std::vector<dso::detail::OceanicTidalWave>::iterator
dso::OceanTide::append_wave(dso::detail::OceanicTidalWave &&wave) {
  if (this->find_tidal_wave(wave.wave().doodson()) != mwaves.end()) {
    /* this doodson number (wave) already exists */
    char buf[8] = "\0";
    fprintf(stderr,
            "[ERROR] Tidal wave with doodson number %s already exists; failed "
            "to append to OceanTide instance (traceback: %s)\n",
            buf, __func__);
    throw std::runtime_error("[ERROR] Failed pushing back tidal wave\n");
  }
  mwaves.emplace_back(wave);
  return mwaves.end() - 1;
}

std::vector<dso::detail::OceanicTidalWave>::iterator
dso::OceanTide::append_wave(const dso::TidalWave &wave, int max_degree,
                            int max_order) {
  if (this->find_tidal_wave(wave.doodson()) != mwaves.end()) {
    /* this doodson number (wave) already exists */
    char buf[8] = "\0";
    fprintf(stderr,
            "[ERROR] Tidal wave with doodson number %s already exists; failed "
            "to append to OceanTide instance (traceback: %s)\n",
            buf, __func__);
    throw std::runtime_error("[ERROR] Failed pushing back tidal wave\n");
  }
  mwaves.emplace_back(
      dso::detail::OceanicTidalWave(wave, max_degree, max_order));
  return mwaves.end() - 1;
}
