#include "earth_rotation.hpp"

int dso::EarthRotation::eops_at(const dso::MjdEpoch &tt) noexcept {
  /* only interpolate if the requested epoch is different than the one we 
   * already have the EOPs for
   */
  if (tt != meop.t()) {
    if (dso::EopSeries::out_of_bounds(meops.interpolate(tt, meop))) {
      fprintf(stderr,
              "[ERROR] Failed interpolating EOPs for given data series "
              "(traceback: %s)\n",
              __func__);
      return 1;
    }
  }
  return 0;
}
