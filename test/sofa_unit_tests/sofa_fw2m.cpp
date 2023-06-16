#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"fw2m"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main() {
  double as[3][3];

  printf("Function #Tests #Fails #Maxerror[sec]    Status  Type\n");
  printf("---------------------------------------------------------------\n");

  int fails = 0;
  double max_error = 0e0;
  for (int i = 0; i < NUM_TESTS; i++) {
    const double gamb = random_angle();
    const double phib = random_angle();
    const double psi = random_angle();
    const double eps = random_angle();
    const auto am = fw2m(gamb, phib, psi, eps);
    iauFw2m(gamb, phib, psi, eps, as);
    if (!approx_equal(am, as)) {
      ++fails;
      /* Angle between rotation matrices [rad] */
      const double theta = rotation_matrix_diff(am, as);
      if (std::abs(theta) > std::abs(max_error)) {
        max_error = theta;
      }
    }
  }
  printf("%8s %6d %6d %+.9e %.7s %s\n", funcs[0], NUM_TESTS, fails,
         dso::rad2sec(max_error), (fails == 0) ? "OK" : "FAILED", "Angle");

  return fails;
}
