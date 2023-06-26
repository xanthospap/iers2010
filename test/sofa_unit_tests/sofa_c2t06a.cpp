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

const char *funcs[] = {"c2t06a"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main() {
  int fails;
  double max_error = 0e0;
  double as[3][3];

  printf("Function #Tests #Fails #Maxerror[masec]  Status  Type\n");
  printf("---------------------------------------------------------------\n");

  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto ut1 = add_random_seconds(tt);
    const double xp = random_angle(-iers2010::DPI / 2, iers2010::DPI / 2);
    const double yp = random_angle(-iers2010::DPI / 2, iers2010::DPI / 2);
    const auto am = c2t06a(tt, ut1, xp, yp);
    iauC2t06a(tt.big() + dso::mjd0_jd, tt.small(), ut1.big() + dso::mjd0_jd,
              ut1.small(), xp, yp, as);
    if (!approx_equal(am, as)) {
      ++fails;
      /* Angle between rotation matrices [rad] */
      const double theta = rotation_matrix_diff(am, as);
      if (std::abs(theta) > std::abs(max_error)) {
        max_error = theta;
      }
    }
  }
  printf("%8s %6d %6d %+.1e %.7s %s\n", funcs[0], NUM_TESTS, fails,
         dso::rad2sec(max_error)*1e3, (fails == 0) ? "OK" : "FAILED",
         "RotMatrix");

  return fails;
}
