#include "unit_test_help.hpp"
#include "iau.hpp"
#include "sofa.h"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 10000;

const char *funcs[] = {"numat"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main(int argc, [[maybe_unused]] char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  int func_it = 0;
  int fails;
  int error = 0;
  double max_error = 0e0;
  double as[3][3];

  printf("Function #Tests #Fails #Maxerror[sec]    Status  Type\n");
  printf("---------------------------------------------------------------\n");

  fails = 0;
  max_error = 0e0;
  for (int i = 0; i < NUM_TESTS; i++) {
    const double epsa = random_angle(-iers2010::DPI / 3, iers2010::DPI / 3);
    const double dpsi = random_angle(-iers2010::DPI / 3, iers2010::DPI / 3);
    const double deps = random_angle(-iers2010::DPI / 3, iers2010::DPI / 3);
    const auto am = numat(epsa, dpsi, deps);
    iauNumat(epsa, dpsi, deps, as);
    if (!approx_equal(am, as)) {
      ++fails;
      /* Angle between rotation matrices [rad] */
      const double theta = rotation_matrix_diff(am, as);
      if (std::abs(theta) > std::abs(max_error)) {
        max_error = theta;
      }
    }
  }
  printf("%8s %6d %6d %+.9e %.7s %s\n", funcs[func_it++], NUM_TESTS, fails,
         dso::rad2sec(max_error), (fails == 0) ? "OK" : "FAILED",
         "RotMatrix");
  if (fails)
    ++error;

  return error;
}
