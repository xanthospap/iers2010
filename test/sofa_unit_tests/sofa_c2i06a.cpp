#include "geodesy/units.hpp"
#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 10000;

const char *funcs[] = {"c2i06a"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main(int argc, [[maybe_unused]] char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  double as[3][3];
  Eigen::Matrix<double,3,3> max_mine, max_sofa;

  printf("Function #Tests #Fails #Maxerror[sec]    Status\n");
  printf("---------------------------------------------------------------\n");

  int fails = 0;
  double max_error = .0e0;
  int i = 0;
  while (i < NUM_TESTS) {
    const auto tt = random_mjd();
    const auto am = num06a(tt);
    iauNum06a(tt.big() + dso::mjd0_jd, tt.small(), as);
    if (!approx_equal(am, as)) {
      ++fails;
      /* angle between rotation matrices [rad] */
      const double theta = rotation_matrix_diff(am, as);
      if (theta == std::numeric_limits<double>::max()) {
        fprintf(stderr, "# Erronuous input @ %s: seti=%d \n", funcs[0], i);
        --i;
      } else {
        if (std::abs(theta) > max_error) {
          max_error = std::abs(theta);
          max_mine = am;
          max_sofa = rxr2mat(as);
        }
      }
    }
    ++i;
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[0], NUM_TESTS, fails,
         dso::rad2sec(max_error), (fails == 0) ? "OK" : "FAILED");

  /*if (fails) {
    // const auto d = max_mine - max_sofa;
    for (int k = 0; k < 3; k++) {
      for (int j = 0; j < 3; j++) {
        printf("%+.12e ", max_mine(k, j));
      }
      printf("\n");
    }
    for (int k = 0; k < 3; k++) {
      for (int j = 0; j < 3; j++) {
        printf("%+.12e ", max_sofa(k, j));
      }
      printf("\n");
    }
  }*/

  return fails;
}
