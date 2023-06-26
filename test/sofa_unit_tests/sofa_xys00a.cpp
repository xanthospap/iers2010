#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 10000;

const char *funcs[] = {"xys00a"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);
const char *args[] = {"xp", "yp", "s"};
const int num_args = sizeof(args) / sizeof(args[0]);

int main() {
  int fails[3] = {0, 0, 0};
  int error = 0;
  double max_error[3] = {0e0, 0e0, 0e0};
  double am[3], as[3];

  printf("Function         #Tests #Fails #Maxerror[masec]  Status  Type\n");
  printf("---------------------------------------------------------------\n");

  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    xys00a(tt, am[0], am[1], am[2]);
    iauXys00a(tt.big() + dso::mjd0_jd, tt.small(), &as[0], &as[1], &as[2]);
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
    printf("%8s %7s %6d %6d %+.1e %.7s %s\n", funcs[0], args[j], NUM_TESTS,
           fails[j], dso::rad2sec(max_error[j])*1e3,
           (fails[j] == 0) ? "OK" : "FAILED", "Angle");
    error += fails[j];
  }

  return error;
}
