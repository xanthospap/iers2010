#include "earth_rotation.hpp"
#include <stdexcept>

dso::EarthRotation::EarthRotation(const char *fn, const dso::MjdEpoch &start_tt,
                                  const dso::MjdEpoch &end_tt) {
  if (dso::parse_iers_C04(fn, start_tt, end_tt, meops)) {
    throw std::runtime_error(
        "[ERROR] Failed constructing EarthRotation instance!\n");
  }
}
