#include "iau.hpp"
#include "unit_test_help.hpp"
#include "geodesy/units.hpp"
#include "sofa.h"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"ee06a", "gmst06", "gmst00", "gst06a", "obl06", "obl80", "ee00", "eect00", "sp00"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main(int argc, char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  int func_it = 0;
  int fails;
  int error = 0;
  double max_error = std::numeric_limits<double>::min();

  printf("Function #Tests #Fails #Maxerror[sec]    Status\n");
  printf("---------------------------------------------------------------\n");

  double am, as;
  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    am = ee06a(tt);
    as = iauEe06a(tt.big()+dso::mjd0_jd, tt.small());
    if (!approx_equal(am, as)) {
      ++fails;
      if (std::abs(am-as) > max_error) max_error = am-as;
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[func_it++], NUM_TESTS, fails, dso::rad2sec(max_error),
         (fails == 0) ? "OK" : "FAILED");
  if (fails) ++error;

  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto ut1 = add_random_seconds(tt);
    am = gmst06(ut1, tt);
    as = iauGmst06(ut1.big() + dso::mjd0_jd, ut1.small(),
                  tt.big() + dso::mjd0_jd, tt.small());
    if (!approx_equal(am, as)) {
      ++fails;
      if (std::abs(am-as) > max_error) max_error = am-as;
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[func_it++], NUM_TESTS, fails, dso::rad2sec(max_error),
         (fails == 0) ? "OK" : "FAILED");
  if (fails) ++error;
  
  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto ut1 = add_random_seconds(tt);
    am = gmst00(ut1, tt);
    as = iauGmst00(ut1.big() + dso::mjd0_jd, ut1.small(),
                  tt.big() + dso::mjd0_jd, tt.small());
    if (!approx_equal(am, as)) {
      ++fails;
      if (std::abs(am-as) > max_error) max_error = am-as;
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[func_it++], NUM_TESTS, fails, dso::rad2sec(max_error),
         (fails == 0) ? "OK" : "FAILED");
  if (fails) ++error;
  
  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto ut1 = add_random_seconds(tt);
    am = gst06a(ut1, tt);
    as = iauGst06a(ut1.big() + dso::mjd0_jd, ut1.small(),
                  tt.big() + dso::mjd0_jd, tt.small());
    if (!approx_equal(am, as)) {
      ++fails;
      if (std::abs(am-as) > max_error) max_error = am-as;
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[func_it++], NUM_TESTS, fails, dso::rad2sec(max_error),
         (fails == 0) ? "OK" : "FAILED");
  if (fails) ++error;
  
  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    am = obl06(tt);
    as = iauObl06(tt.big()+dso::mjd0_jd, tt.small());
    if (!approx_equal(am, as)) {
      ++fails;
      if (std::abs(am-as) > max_error) max_error = am-as;
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[func_it++], NUM_TESTS, fails, dso::rad2sec(max_error),
         (fails == 0) ? "OK" : "FAILED");
  if (fails) ++error;
  
  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    am = obl80(tt);
    as = iauObl80(tt.big()+dso::mjd0_jd, tt.small());
    if (!approx_equal(am, as)) {
      ++fails;
      if (std::abs(am-as) > max_error) max_error = am-as;
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[func_it++], NUM_TESTS, fails, dso::rad2sec(max_error),
         (fails == 0) ? "OK" : "FAILED");
  if (fails) ++error;
  
  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const double epsa = random_angle(-iers2010::DPI/3, iers2010::DPI/3);
    const double dpsi = random_angle(-iers2010::DPI/3, iers2010::DPI/3);
    am = ee00(tt, epsa, dpsi);
    as = iauEe00(tt.big()+dso::mjd0_jd, tt.small(), epsa, dpsi);
    if (!approx_equal(am, as)) {
      ++fails;
      if (std::abs(am-as) > max_error) max_error = am-as;
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[func_it++], NUM_TESTS, fails, dso::rad2sec(max_error),
         (fails == 0) ? "OK" : "FAILED");
  if (fails) ++error;
  
  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    am = eect00(tt);
    as = iauEect00(tt.big()+dso::mjd0_jd, tt.small());
    if (!approx_equal(am, as)) {
      ++fails;
      if (std::abs(am-as) > max_error) max_error = am-as;
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[func_it++], NUM_TESTS, fails, dso::rad2sec(max_error),
         (fails == 0) ? "OK" : "FAILED");
  if (fails) ++error;
  
  fails = 0;
  max_error = std::numeric_limits<double>::min();
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    am = sp00(tt);
    as = iauSp00(tt.big()+dso::mjd0_jd, tt.small());
    if (!approx_equal(am, as)) {
      ++fails;
      if (std::abs(am-as) > max_error) max_error = am-as;
    }
  }
  printf("%8s %6d %6d %+.9e %s\n", funcs[func_it++], NUM_TESTS, fails, dso::rad2sec(max_error),
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
