#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"ee06a"};
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

  double am, as;
  fails = 0;
  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    // add a few seconds to TT to get a random UT1 date
    auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();

    am = ee06a(tt);
    as = iauEe06a(jdtt._big, jdtt._small);
    if (!approx_equal(am, as))
      ++fails;
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
