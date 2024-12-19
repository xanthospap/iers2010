#include "atmospheric_tides.hpp"
#include <array>
#include <stdexcept>

namespace {} /* anonymous namespace */

const dso::detail::TidalConstituentArrayEntry *
dso::find_wave_entry(const char *name) noexcept {
  const auto it =
      std::find_if(dso::detail::AtmosphericTidalHarmonics.begin(),
                   dso::detail::AtmosphericTidalHarmonics.end(),
                   [=](const dso::detail::TidalConstituentArrayEntry &wave) {
                     return !std::strcmp(wave.name(), name);
                   });
  return (it == dso::detail::AtmosphericTidalHarmonics.end()) ? nullptr
                                                              : &(*it);
}

dso::TidalWave dso::get_wave(const char *name) {
  const auto it = dso::find_wave_entry(name);
  if (it)
    return dso::TidalWave(*it);
  throw std::runtime_error(
      "[ERROR] Failed matching tidal wave name to AtmosphericTidalHarmonics\n");
}
