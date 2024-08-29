#include "doodson.hpp"
#include <cstring>

const dso::detail::TidalConstituentsArrayEntry *
dso::find_wave_entry(const char *name) noexcept {
  for (const auto &w : dso::AtmosphericTidalHarmonics) {
    if (!std::strcmp(name, w._n))
      return &w;
  }
  // throw std::runtime_error("[ERROR] failed finding Doodson wave!\n");
  return nullptr;
}

dso::TidalWave
dso::get_wave(const char *name) noexcept {
  auto w = dso::find_wave_entry(name);
  if (w) return dso::TidalWave(w);
  throw std::runtime_error("[ERROR] failed finding Doodson wave!\n");
}
