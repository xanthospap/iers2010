#include "iers2010.hpp"
#include <cassert>

int main() {
  printf("Testing Lagrangian interpolation against pre-computed results\n");

  constexpr const int N = 10;
  const double x0[N] = {
    0.5e0, 0.6e0, 0.99e0, 1.2e0, 2e0, 3e0, 4e0, 5e0, 6e0, 7e0};
  const double y0[N] = {
    99e0, 99e0, 99e0, 99e0, 2e0, 6e0, 24e0, 120e0, 99e0, 99e0};

  double yout;
  int k = -1;
  int error;
  error = iers2010::interp::lagint(x0, y0, N, 3.289e0, yout, k, 3);
  assert(!error);
  printf("X=%.3f Y=%.3f idx=%d\n", 3.289, yout, k);

  return 0;
}
