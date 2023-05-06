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

const char *funcs[] = {"nut00a"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);
const char *args[] = {"dpsi", "deps"};
const int num_args = sizeof(args) / sizeof(args[0]);

int main(int argc, [[maybe_unused]] char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  int fails[2]={0,0};
  int error = 0;
  double max_error[2] = {0e0,0e0};
  double am[2], as[2];

  printf("Function         #Tests #Fails #Maxerror[sec]    Status\n");
  printf("---------------------------------------------------------------\n");

  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    nut00a(tt, am[0], am[1]);
    iauNut00a(tt.big() + dso::mjd0_jd, tt.small(), &as[0], &as[1]);
    for (int j = 0; j < 2; j++) {
      if (!approx_equal(am[j], as[j])) {
        ++fails[j];
        if (std::abs(am[j] - as[j]) > max_error[j]) {
          max_error[j] = std::abs(am[j] - as[j]);
        }
      }
    }
  }

  for (int j = 0; j < 2; j++) {
    printf("%8s %7s %6d %6d %+.9e %s\n", funcs[0], args[j], NUM_TESTS,
           fails[j], dso::rad2sec(max_error[j]), (fails[j] == 0) ? "OK" : "FAILED");
  }

  error=0;
  for (int j=0;j<2;j++) error+=fails[j];

  return error;
}
