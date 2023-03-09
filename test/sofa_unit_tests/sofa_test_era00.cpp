#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <datetime/dtcalendar.hpp>
#include <datetime/dtfund.hpp>
#include <geodesy/units.hpp>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 1000000;

const char *funcs[] = {"era00"};
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
  fprintf(stderr, "%8s %12s\n", "Era [deg]", "dEra [sec]");
  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    am = era00(tt);

    auto foo = tt.jd_split<dso::TwoPartDate::JdSplitMethod::DT>();
    as = iauEra00(foo._big, foo._small);
    if (!approx_equal(am, as)){
      fprintf(stderr, "%.12e %+9.4f %.12e\n", tt.mjd(), dso::rad2deg(am),
              dso::rad2sec(am - as));
      ++fails;
    }
  }

  // how much of an angle is 1sec?
  {
    const int tests = 10000;
    double fac[] = {1e0, 1e-3, 1e-6, 1e-9};
    for (int k = 0; k < 4; k++) {
      double sum = 0e0;
      const double oneSec = fac[k] / 86400e0;
      for (int i = 0; i < tests; i++) {
        const auto tt = random_mjd();
        am = era00(tt);
        const dso::TwoPartDate ttp(tt._big, tt._small + oneSec);
        as = era00(ttp);
        sum += dso::rad2sec(as - am);
      }
      printf("%.2esec is about ~%.12e [arcsec] ERA angle\n", fac[k],
             sum / tests);
    }
  }

  /*
  double jd_start = 2460004.5e0;
  double cjd = jd_start;
  double seconds = 0e0;
  while (cjd<jd_start + 2.1) {
    const double fday = seconds / 86400e0;
    
    const auto mjd = dso::TwoPartDate(cjd-dso::mjd0_jd, fday);
    am = era00(mjd);

    const auto jd = dso::TwoPartDate(cjd, fday);
    as = iauEra00(jd._big, jd._small);
    if (!approx_equal(am, as)){
      fprintf(stderr, "%+9.4f %.12e\n", dso::rad2deg(am), dso::rad2sec(am-as));
      ++fails;
    }

    seconds += 1e0;
    if (seconds > 86400e0){
      seconds = 0;
      cjd += 1e0;
    }
  }
  */
  
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
