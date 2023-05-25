#include "geodesy/units.hpp"
#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <datetime/dtcalendar.hpp>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"c2ixys"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main(int argc, [[maybe_unused]] char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  int func_it = 0;
  double as[3][3];

  printf("Function #Tests #Fails #Maxerror[sec]    Status\n");
  printf("---------------------------------------------------------------\n");

  int fails = 0;
  double max_error = 0e0;
  for (int i = 0; i < NUM_TESTS; i++) {
    const double xp = random_angle(-iers2010::DPI / 2, iers2010::DPI / 2);
    const double yp = random_angle(-iers2010::DPI / 2, iers2010::DPI / 2);
    const double s = random_angle(-iers2010::DPI / 2, iers2010::DPI / 2);
    const auto am = c2ixys(xp, yp, s);
    iauC2ixys(xp, yp, s, as);
    if (!approx_equal(am, as)) {
      ++fails;
      /* angle between rotation matrices [rad] */
      const double theta = rotation_matrix_diff(am, as);
      if (std::abs(theta) > max_error) {
        max_error = theta;
      }
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[func_it++], NUM_TESTS, fails,
         dso::rad2sec(max_error), (fails == 0) ? "OK" : "FAILED");

  return fails;
}
