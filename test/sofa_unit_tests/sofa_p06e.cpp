#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"p06e"};
const char *args[] = {"eps0", "psia", "oma",  "bpa", "bqa",   "pia",
                      "bpia", "epsa", "chia", "za",  "zetaa", "thetaa",
                      "pa",   "gam",  "phi",  "psi"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);
const int num_args = sizeof(args) / sizeof(args[0]);

int main() {
  int fails[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int error = 0;
  double max_error[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double as[16], am[16];

  printf("Function       #Tests #Fails #Maxerror[sec]    Status  Type\n");
  printf("---------------------------------------------------------------\n");

  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    p06e(tt, am[0], am[1], am[2], am[3], am[4], am[5], am[6], am[7], am[8],
         am[9], am[10], am[11], am[12], am[13], am[14], am[15]);
    iauP06e(tt.big() + dso::mjd0_jd, tt.small(), &as[0], &as[1], &as[2], &as[3],
            &as[4], &as[5], &as[6], &as[7], &as[8], &as[9], &as[10], &as[11],
            &as[12], &as[13], &as[14], &as[15]);
    for (int j = 0; j < 16; j++) {
      if (!approx_equal(am[j], as[j])) {
        ++fails[j];
        /* Angle between rotation matrices [rad] */
        const double diff = as[j] - am[j];
        if (std::abs(diff) > std::abs(max_error[j])) {
          max_error[j] = diff;
        }
      }
    }
  }

  for (int i = 0; i < 16; i++) {
    printf("%8s %8s %6d %6d %+.9e %.7s %s\n", funcs[0], args[i], NUM_TESTS,
           fails[i], dso::rad2sec(max_error[i]),
           (fails[i] == 0) ? "OK" : "FAILED", "Angle");
  }

  error = 0;
  for (int i = 0; i < 16; i++)
    error += fails[i];

  return error;
}
