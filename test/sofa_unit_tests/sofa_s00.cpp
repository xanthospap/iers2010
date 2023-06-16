#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"s00"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main() {
  int fails;
  double max_error = 0e0;
  double am, as;

  printf("Function         #Tests #Fails #Maxerror[sec]    Status  Type\n");
  printf("---------------------------------------------------------------\n");

  fails = 0;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const double xp = random_angle(-iers2010::DPI / 4, iers2010::DPI / 4);
    const double yp = random_angle(-iers2010::DPI / 4, iers2010::DPI / 4);
    am = s00(tt, xp, yp);
    as = iauS00(tt.big() + dso::mjd0_jd, tt.small(), xp, yp);
    if (!approx_equal(am, as)) {
      ++fails;
      if (std::abs(am - as) > std::abs(max_error)) {
        max_error = as - am;
      }
    }
  }

  printf("%8s %6d %6d %+.9e %.7s %s\n", funcs[0], NUM_TESTS, fails,
         dso::rad2sec(max_error), (fails == 0) ? "OK" : "FAILED", "Angle");

  return fails;
}
