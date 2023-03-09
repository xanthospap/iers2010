#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 1000;

const char *funcs[] = {"xy06"};
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
  fprintf(stderr, "%7s %7s       %12s %12s\n", "X[deg]", "Y[deg]", "dX[asec]",
          "dY[asec]");

  double am[2], as[2];
  fails = 0;
  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    xy06(tt, am[0], am[1]);

    // transform the MJD to a JD split using the J2000 method
    auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();
    assert(jdtt._big == dso::j2000_jd);
    iauXy06(jdtt._big, jdtt._small, as, as+1);
    if (!approx_equal(am[0], as[0]) || !approx_equal(am[1], as[1])) {
      ++fails;
      fprintf(stderr, "%+7.3f %+7.3f [deg] %+.6e %+.6e [asec]\n",
              dso::rad2deg(am[0]), dso::rad2deg(am[1]),
              dso::rad2sec(am[0] - as[0]), dso::rad2sec(am[1] - as[1]));
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
