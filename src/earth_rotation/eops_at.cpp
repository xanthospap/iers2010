#include "earth_rotation.hpp"

int dso::EarthRotation::eops_at(const dso::MjdEpoch &tt) noexcept {
  if (tt != meop.t()) {
    const auto res = meops.interpolate(tt, meop);
    switch (res) {
    case dso::EopSeries::EopInterpolationResult::Linear:
    case dso::EopSeries::EopInterpolationResult::PolyDegreeDescreased:
    case dso::EopSeries::EopInterpolationResult::PolyDegreeRequested:
      return 0;
    case dso::EopSeries::EopInterpolationResult::OutOfBoundsPrior:
    case dso::EopSeries::EopInterpolationResult::OutOfBoundsLater:
      fprintf(stderr,
              "[ERROR] Failed interpolating EOPs for given data series "
              "(traceback: %s)\n",
              __func__);
      return 1;
    }
  }
  return 0;
}
