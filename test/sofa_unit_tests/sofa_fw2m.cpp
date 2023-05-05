#include "iau.hpp"
#include "unit_test_help.hpp"
#include "geodesy/units.hpp"
#include "sofa.h"
#include <cassert>
#include <cstdio>
#include <datetime/dtcalendar.hpp>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 10000;

const char *funcs[] = {"fw2m"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main(int argc, [[maybe_unused]]char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  int func_it = 0;
  int fails;
  int error = 0;
  double max_error = std::numeric_limits<double>::min();
  double as[3][3];

  printf("Function #Tests #Fails #Maxerror[sec]    Status\n");
  printf("---------------------------------------------------------------\n");

  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const double gamb = random_angle();
    const double phib = random_angle();
    const double psi = random_angle();
    const double eps = random_angle();
    const auto am = fw2m(gamb, phib, psi, eps);
    iauFw2m(gamb, phib, psi, eps, as);
    if (!approx_equal(am, as)) {
      ++fails;
      /* angle between rotation matrices [rad] */
      const double theta = rotation_matrix_diff(am, as);
      if (std::abs(theta) > max_error) {
        max_error = theta;
      }
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[func_it++], NUM_TESTS, fails, dso::rad2sec(max_error),
         (fails == 0) ? "OK" : "FAILED");
  if (fails) ++error;

  return error;
}
