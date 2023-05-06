#include "iau.hpp"
#include "unit_test_help.hpp"
#include "geodesy/units.hpp"
#include "sofa.h"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"era00"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main() {
  int fails;
  int error = 0;
  double max_error = std::numeric_limits<double>::min();

  printf("Function #Tests #Fails #Maxerror[sec]    Status\n");
  printf("---------------------------------------------------------------\n");

  double am, as;
  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto ut = random_mjd();
    am = era00(ut);
    as = iauEra00(ut.big()+dso::mjd0_jd, ut.small());
    if (!approx_equal(am, as)) {
      ++fails;
      if (std::abs(am-as) > max_error) max_error = am-as;
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[0], NUM_TESTS, fails, dso::rad2sec(max_error),
         (fails == 0) ? "OK" : "FAILED");
  if (fails) ++error;
  
  return error;
}
