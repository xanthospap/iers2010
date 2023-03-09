#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"s06"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main(int argc, char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  int func_it = 0;
  int fails;
  int error = 0;

  printf("Function #Tests #Fails Status\n");
  printf("---------------------------------------------------------------\n");
  fprintf(stderr, "%7s %12s\n", "S[deg]", "dS[asec]");

  double am, as;
  fails = 0;
  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    // random CIP components, small angles in X, Y [rad]
    const double X = random_angle(-dso::deg2rad(.9e0), dso::deg2rad(.9e0));
    const double Y = random_angle(-dso::deg2rad(.09e0), dso::deg2rad(.09e0));
    am = s06(tt, X, Y);

    // transform the MJD to a JD split using the J2000 method
    auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();
    assert(jdtt._big == dso::j2000_jd);
    as = iauS06(jdtt._big, jdtt._small, X, Y);
    if (!approx_equal(am, as)) {
      ++fails;
      fprintf(stderr, "%+7.3f %+.6e\n", dso::rad2deg(am),
              dso::rad2sec(am - as));
    }
  }
  printf("%8s %6d %6d %s\n", funcs[func_it++], NUM_TESTS, fails,
         (fails == 0) ? "OK" : "FAILED");
  if (fails) ++error;

  assert(num_funs == func_it);

  if (!error)
    printf("Program %s : All tests passed!\n", argv[0]);
  else
    printf("Program %s : CHECK FAILED!\n", argv[0]);
  printf("---------------------------------------------------------------\n");
  return error;
}
