#include "doodson.hpp"

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#include <strings.h>
#elif defined _MSC_VER
#define strcasecmp _stricmp
#endif

const dso::detail::TidalConstituentArrayEntry *
match_ocean_tide_wave(const char *wave) noexcept {
  for (const auto &w : dso::detail::OceanicTidalHarmonics) {
    const char *a = w.name();
    const char *b = wave;
    if (!strcasecmp(a, b))
      return &w;
  }
  return nullptr;
}

const dso::detail::TidalConstituentArrayEntry *
match_ocean_tide_wave(const dso::DoodsonConstituent &wave) noexcept {
  for (const auto &w : dso::detail::OceanicTidalHarmonics) {
    if (wave == w._d) return &w;
  }
  return nullptr;
}
