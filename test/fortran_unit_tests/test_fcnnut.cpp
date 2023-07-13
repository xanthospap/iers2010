#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>
#include <cstdio>
#include <datetime/dtcalendar.hpp>
#include <limits>

using namespace iers2010;

constexpr const int NUM_TESTS = 10000;

const char *funcs[] = {"fcnnut"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main() {

  double x, y, dx, dy;
  for (int i = 0; i < NUM_TESTS; i++) {
    /* random MJD in range 1980 to 2030 */
    double mjd_double = generate_random_double(44239e0, 62502e0);
    dso::TwoPartDate mjd_tt(mjd_double, 0e0);
    fcnnut(mjd_tt, x, y, dx, dy);
    printf("%+.20e %+.20e %+.20e %+.20e %+.20e\n", mjd_double, x, y, dx, dy);
  }

  return 0;
}
