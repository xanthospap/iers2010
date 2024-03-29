#include "iers2010.hpp"

namespace {
constexpr double orthow[3][12] = {
    {-6.77832e0, -14.86323e0, 0.47884e0, -1.45303e0, 0.16406e0, 0.42030e0,
     0.09398e0, 25.73054e0, -4.77974e0, 0.28080e0, 1.94539e0, -0.73089e0},
    {14.86283e0, -6.77846e0, 1.45234e0, 0.47888e0, -0.42056e0, 0.16469e0,
     15.30276e0, -4.30615e0, 0.07564e0, 2.28321e0, -0.45717e0, -1.62010e0},
    {-1.76335e0, 1.03364e0, -0.27553e0, 0.34569e0, -0.12343e0, -0.10146e0,
     -0.47119e0, 1.28997e0, -0.19336e0, 0.02724e0, 0.08955e0, 0.04726e0}};
}

int iers2010::ortho_eop(const dso::TwoPartDate &t, double &dx, double &dy,
                        double &dut1) noexcept {

  // Compute the partials of the tidal variations to the orthoweights
  double h[12];
  iers2010::oeop::cnmtx(t, h);

  // Compute eop changes
  dx = dy = dut1 = 0e0;

  for (int j = 0; j < 12; j++) {
    dx += h[j] * orthow[0][j];
    dy += h[j] * orthow[1][j];
    dut1 += h[j] * orthow[2][j];
  }

  // finished
  return 0;
}
