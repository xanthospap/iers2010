#include "tide_atlas.hpp"

std::vector<dso::TidalConstituent>::iterator
dso::TideAtlas::append_wave(const dso::TidalConstituent &wave) {
  if (this->find_tidal_wave(wave.wave().doodson()) != mwaves.end()) {
    /* this doodson number (wave) already exists */
    char buf[8] = "\0";
    fprintf(stderr,
            "[ERROR] Tidal wave with doodson number %s already exists; failed "
            "to append to TideAtlas instance (traceback: %s)\n",
            wave.wave().doodson().str(buf, true), __func__);
    throw std::runtime_error("[ERROR] Failed pushing back tidal wave\n");
  }
  mwaves.emplace_back(wave);
  return mwaves.end() - 1;
}

std::vector<dso::TidalConstituent>::iterator
dso::TideAtlas::append_wave(const dso::TidalWave &wave, int max_degree,
                            int max_order) {
  if (this->find_tidal_wave(wave.doodson()) != mwaves.end()) {
    /* this doodson number (wave) already exists */
    char buf[8] = "\0";
    fprintf(stderr,
            "[ERROR] Tidal wave with doodson number %s already exists; failed "
            "to append to TideAtlas instance (traceback: %s)\n",
            wave.doodson().str(buf, true), __func__);
    throw std::runtime_error("[ERROR] Failed pushing back tidal wave\n");
  }
  mwaves.emplace_back(
      dso::TidalConstituent(wave, max_degree, max_order));
  return mwaves.end() - 1;
}
