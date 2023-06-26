#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"pom00"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main() {
  int error = 0;
  double as[3][3];

  printf("Function #Tests #Fails #Maxerror[masec]  Status  Type\n");
  printf("---------------------------------------------------------------\n");

  int fails = 0;
  double max_error = 0e0;
  for (int i = 0; i < NUM_TESTS; i++) {
    const double xp = random_angle(-iers2010::DPI / 4, iers2010::DPI / 4);
    const double yp = random_angle(-iers2010::DPI / 4, iers2010::DPI / 4);
    const double s = random_angle(-iers2010::DPI / 3, iers2010::DPI / 3);
    const auto am = pom00(xp, yp, s);
    iauPom00(xp, yp, s, as);
    if (!approx_equal(am, as)) {
      ++fails;
      /* Angle between rotation matrices [rad] */
      const double theta = rotation_matrix_diff(am, as);
      if (std::abs(theta) > std::abs(max_error)) {
        max_error = theta;
      }
    }
  }
  printf("%8s %6d %6d %+.1e %.7s %s\n", funcs[0], NUM_TESTS, fails,
         dso::rad2sec(max_error)*1e3, (fails == 0) ? "OK" : "FAILED",
         "RotMatrix");
  if (fails)
    ++error;

  return error;
}
