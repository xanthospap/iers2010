#include "tropo.hpp"

double dso::asknewet(double e, double Tm, double lambda) noexcept {
  // coefficients
  constexpr double k1 = 77.604e0;                         // K / hPa
  constexpr double k2 = 64.79e0;                          // K / hPa
  constexpr double k2p = k2 - k1 * 18.0152e0 / 28.9644e0; // K / hPa
  constexpr double k3 = 377600e0;                         // KK / hPa

  // mean gravity in m / s ** 2
  constexpr double gm = 9.80665e0;
  // molar mass of dry air in kg / mol
  constexpr double dMtr = 28.9651e-3;
  // universal gas constant in J / K / mol
  constexpr double R = 8.3143e0;

  // specific gas constant for dry consituents
  constexpr double Rd = R / dMtr;

  return 1e-6 * (k2p + k3 / Tm) * Rd / (lambda + 1) / gm * e;
}
