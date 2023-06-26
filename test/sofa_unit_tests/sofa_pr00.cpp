#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 10000;

const char *funcs[] = {"pr00"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);
const char *args[] = {"dpsipr", "depspr"};
const int num_args = sizeof(args) / sizeof(args[0]);

int main(int argc, [[maybe_unused]] char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  int fails[] = {0, 0};
  int error = 0;
  double max_error[] = {0e0, 0e0};
  double sdpsi, sdeps, mdpsi, mdeps;

  printf("Function         #Tests #Fails #Maxerror[masec]  Status  Type\n");
  printf("---------------------------------------------------------------\n");

  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    pr00(tt, mdpsi, mdeps);
    iauPr00(tt.big() + dso::mjd0_jd, tt.small(), &sdpsi, &sdeps);
    /* dpsipr */
    if (!approx_equal(sdpsi, mdpsi)) {
      ++fails[0];
      if (std::abs(sdpsi - mdpsi) > std::abs(max_error[0])) {
        max_error[0] = sdpsi - mdpsi;
      }
    }
    /* deps */
    if (!approx_equal(sdeps, mdeps)) {
      ++fails[1];
      if (std::abs(sdeps - mdeps) > std::abs(max_error[1])) {
        max_error[1] = sdeps - mdeps;
      }
    }
  }

  for (int j = 0; j < 2; j++) {
    printf("%8s %7s %6d %6d %+.1e %.7s %s\n", funcs[0], args[j], NUM_TESTS,
           fails[j], dso::rad2sec(max_error[j])*1e3,
           (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  }

  error = 0;
  for (int j = 0; j < 2; j++)
    error += fails[j];

  return error;
}
