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

const char *funcs[] = {"s06"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main() {
  int fails;
  int error = 0;
  double max_error = 0e0;
  double am, as;

  printf("Function         #Tests #Fails #Maxerror[sec]    Status\n");
  printf("---------------------------------------------------------------\n");

  fails = 0;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const double xp = random_angle(-iers2010::DPI/2, iers2010::DPI/2);
    const double yp = random_angle(-iers2010::DPI/2, iers2010::DPI/2);
    am = s06(tt, xp, yp);
    as = iauS06(tt.big() + dso::mjd0_jd, tt.small(), xp, yp);
      if (!approx_equal(am, as)) {
        ++fails;
        if (std::abs(am - as) > max_error) {
          max_error = std::abs(am - as);
        }
      }
    }

    printf("%8s %6d %6d %+.9e %s\n", funcs[0], NUM_TESTS,
           fails, dso::rad2sec(max_error), (fails == 0) ? "OK" : "FAILED");
  if (fails)
    ++error;

  return error;
}
