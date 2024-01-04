#include "doodson.hpp"
#include <stdexcept>
#include <cstring>

dso::detail::TidalConstituentsArrayEntry dso::get_wave(const char *name) {
  for (const auto &w : dso::TidalConstituentsArray) {
    if (!std::strcmp(name, w._n))
      return w;
  }
  throw std::runtime_error("[ERROR] failed finding Doodson wave!\n");
}
