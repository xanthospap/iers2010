#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"era00"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main() {

  printf("Function #Tests #Fails #Maxerror[sec]    Status  Type\n");
  printf("---------------------------------------------------------------\n");

  //double mean_l=0e0, mean_r=0e0;
  //double max_l=0e0, max_r=0e0;

  double am, as;
  int fails = 0;
  double max_error = 0e0;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto ut = random_mjd();
    am = era00(ut);
    as = iauEra00(ut.big() + dso::mjd0_jd, ut.small());
    if (!approx_equal(am, as)) {
      ++fails;
      // printf("diff = %.20f - %.20f = %.20f tt=%.3f + %.20f\n", as, am, as-am, ut.big(), ut.small());
      if (std::abs(am - as) > std::abs(max_error)) {
        max_error = as - am;
      }
      //if (ut.small()<.5e0) {
      //  mean_l += std::abs(am - as);
      //  if (std::abs(as-am)>max_l) max_l = as-am;
      //} else {
      //  mean_r += std::abs(am-as);
      //  if (std::abs(as-am)>max_r) max_r = as-am;
      //}
    }
  }
  printf("%8s %6d %6d %+.9e %.7s %s\n", funcs[0], NUM_TESTS, fails,
         dso::rad2sec(max_error), (fails == 0) ? "OK" : "FAILED", "Angle");

  //printf("Means: left=%.9e right=%.9e\n", dso::rad2sec(mean_l) / NUM_TESTS,
  //       dso::rad2sec(mean_r) / NUM_TESTS);
  //printf("Max  : left=%.9e right=%.9e\n", dso::rad2sec(max_l),
  //       dso::rad2sec(max_r));
  return fails;
}
