#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 1000000;

const char *funcs[] = {"c2ixys"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main(int argc, [[maybe_unused]] char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  double as[3][3];

  printf("Function #Tests #Fails #Maxerror[sec]    Status\n");
  printf("---------------------------------------------------------------\n");

  int fails = 0;
  double max_error = 0e0;
  int i = 0;
  while (i < NUM_TESTS) {
    const double xp = random_angle(-iers2010::DPI / 5, iers2010::DPI / 5);
    const double yp = random_angle(-iers2010::DPI / 5, iers2010::DPI / 5);
    const double s = random_angle(-iers2010::DPI / 5, iers2010::DPI / 5);
    /* check for erronuous input */
    int error_in = 0;
    {
      const double r2 = xp * xp + yp * yp;
      const double e = (r2 > 0e0) ? std::atan2(yp, xp) : 0e0;
      const double d = std::atan(std::sqrt(r2 / (1e0 - r2)));
      if (e != e || d != d)
        ++error_in;
    }
    /* if input is ok */
    if (!error_in) {
      const auto am = c2ixys(xp, yp, s);
      iauC2ixys(xp, yp, s, as);
      if (!approx_equal(am, as)) {
        ++fails;
        /* Angle between rotation matrices [rad] */
        const double theta = rotation_matrix_diff(am, as);
        if (theta == std::numeric_limits<double>::max()) {
          --i;
        } else {
          if (std::abs(theta) > std::abs(max_error)) {
            max_error = theta;
          }
        }
      }
      ++i;
    } else {
      /* erronuous input */
      fprintf(stderr, "# Erronuous input @ %s: (%.12e %.12e %.12e) seti=%d \n",
              funcs[0], xp, yp, s, i);
    }
  }
  printf("%8s %6d %6d %+.9e %.7s %s\n", funcs[0], NUM_TESTS, fails,
         dso::rad2sec(max_error), (fails == 0) ? "OK" : "FAILED",
         "RotMatrix");

  return fails;
}
