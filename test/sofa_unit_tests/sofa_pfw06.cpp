#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 10000;

const char *funcs[] = {"pfw06"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);
const char *args[] = {"gamb", "phib", "psib", "epsa"};
const int num_args = sizeof(args) / sizeof(args[0]);

int main() {
  int fails[4] = {0, 0, 0, 0};
  int error = 0;
  double max_error[4] = {0, 0, 0, 0};
  double am[4], as[4];

  printf("Function         #Tests #Fails #Maxerror[masec]  Status  Type\n");
  printf("---------------------------------------------------------------\n");

  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    pfw06(tt, am[0], am[1], am[2], am[3]);
    iauPfw06(tt.big() + dso::mjd0_jd, tt.small(), &as[0], &as[1], &as[2],
             &as[3]);
    for (int j = 0; j < 4; j++) {
      if (!approx_equal(am[j], as[j])) {
        ++fails[j];
        if (std::abs(am[j] - as[j]) > std::abs(max_error[j])) {
          max_error[j] = as[j] - am[j];
        }
      }
    }
  }

  for (int j = 0; j < 4; j++) {
    printf("%8s %7s %6d %6d %+.1e %.7s %s\n", funcs[0], args[j], NUM_TESTS,
           fails[j], dso::rad2sec(max_error[j])*1e3,
           (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  }

  error = 0;
  for (int j = 0; j < 4; j++)
    error += fails[j];

  return error;
}
