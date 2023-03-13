#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>

/*
 * check conversion of a random date (assuming it is Julian), to Julian 
 * Centuries since J2000
 * This function is heavily used in SOFA and this library and results should
 * match
 */

constexpr const double DJC (36525e0);
constexpr const double DJM (365250e0);
constexpr const double DJ00 (2451545e0);

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"-"};
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
    // my result for julian centuries since J2000
    am = tt.jcenturies_sinceJ2000();
    
    // how i would split a JD, using the J2000 method (SOFA)
    auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();
    assert(jdtt._big == dso::j2000_jd);
    // ... and then SOFA does this computation to compute JC since J2000
    as = ((jdtt._big - DJ00) + jdtt._small) / DJC;

    // the two Julian Centuries should be the same ...
    if (!approx_equal(am, as))
      ++fails;
    // printf("%.15e %.15e\n", jdtt.mjd(), am-as);
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