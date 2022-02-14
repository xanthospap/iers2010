#include "tropo.hpp"
#include <cmath>

double dso::saasthyd(double p, double dlat, double hell) noexcept {
  // calculate denominator f
  const double f = 1e0 - 0.00266e0 * std::cos(2e0 * dlat) - 0.00000028e0 * hell;

  // calculate the zenith hydrostatic delay zhd
  return 0.0022768 * p / f;
}
