#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 10000;

const char *funcs[] = {"xys06a"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);
const char *args[] = {"xp", "yp", "s"};
const int num_args = sizeof(args) / sizeof(args[0]);

int main() {
  int fails[3] = {0, 0, 0};
  int error = 0;
  double max_error[3] = {0, 0, 0};
  double am[3], as[3];

  printf("Function         #Tests #Fails #Maxerror[sec]    Status\n");
  printf("---------------------------------------------------------------\n");

  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    xys06a(tt, am[0], am[1], am[2]);
    iauXys06a(tt.big() + dso::mjd0_jd, tt.small(), &as[0], &as[1], &as[2]);
    for (int j = 0; j < 3; j++) {
      if (!approx_equal(am[j], as[j])) {
        ++fails[j];
        if (std::abs(am[j] - as[j]) > std::abs(max_error[j])) {
          max_error[j] = as[j] - am[j];
        }
      }
    }
  }

  error = 0;
  for (int j = 0; j < 3; j++) {
    printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[0], args[j], NUM_TESTS,
           fails[j], dso::rad2sec(max_error[j]),
           (fails[j] == 0) ? "OK" : "FAILED", "Angle");
    error += fails[j];
  }

  return error;
}
