#include "geodesy/units.hpp"
#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <datetime/dtcalendar.hpp>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 10000;

const char *funcs[] = {"xy06"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);
const char *args[] = {"xp", "yp"};
const int num_args = sizeof(args) / sizeof(args[0]);

int main(int argc, [[maybe_unused]] char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  int func_it = 0;
  int fails;
  int error = 0;
  double max_error[2] = {std::numeric_limits<double>::min()};
  double am[2], as[2];

  printf("Function         #Tests #Fails #Maxerror[sec]    Status\n");
  printf("---------------------------------------------------------------\n");

  fails = 0;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    xy06(tt, am[0], am[1]);
    iauXy06(tt.big() + dso::mjd0_jd, tt.small(), &as[0], &as[1]);
    for (int j = 0; j < 2; j++) {
      if (!approx_equal(am[j], as[j])) {
        ++fails;
        if (std::abs(am[j] - as[j]) > max_error[j]) {
          max_error[j] = std::abs(am[j] - as[j]);
        }
      }
    }
  }

  for (int j = 0; j < 2; j++) {
    printf("%8s/%7s %6d %6d %+.9e %s\n", funcs[func_it], args[j], NUM_TESTS,
           fails, dso::rad2sec(max_error[j]), (fails == 0) ? "OK" : "FAILED");
  }
  if (fails)
    ++error;

  return error;
}
