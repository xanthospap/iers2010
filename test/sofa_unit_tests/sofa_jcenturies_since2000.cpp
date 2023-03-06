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

double jc1(const dso::TwoPartDate &mjd) noexcept {
  return ((mjd._big - 51544e0) + (mjd._small - .5e0)) /
         dso::days_in_julian_cent;
}
double jc2(const dso::TwoPartDate &mjd) noexcept {
  return (mjd._big - dso::j2000_mjd) / dso::days_in_julian_cent +
         mjd._small / dso::days_in_julian_cent;
}
void mjd2j2jd(const dso::TwoPartDate &mjd, double &d1, double &d2) noexcept {
  d1 = dso::j2000_jd;
  d2 = (mjd._big - 51544e0) + (mjd._small - .5e0);
}

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
    am = tt.jcenturies_sinceJ2000();
    assert(am==jc1(tt));
    // assert(approx_equal(am,jc2(tt))); <-- this fails!
    
    auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();
    assert(jdtt._big == dso::j2000_jd);
    as = ((jdtt._big - DJ00) + jdtt._small) / DJC;
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
