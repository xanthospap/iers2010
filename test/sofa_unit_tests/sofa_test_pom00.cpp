#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"pom00"};
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
  fprintf(stderr, "Unequal Components: #/9 Max Diff (Abs Value) Value at Max Dif\n");

  fails = 0;
  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    // get s' for date [rad]
    const double sp = sp00(tt);
    // random values for xp, yp [-2, 2] degrees
    const double xp = random_angle(-dso::deg2rad(2), dso::deg2rad(2));
    const double yp = random_angle(-dso::deg2rad(2), dso::deg2rad(2));

    // get matrix
    const auto am = pom00(xp,yp,sp);

    // get SOFA matrix
    double as[3][3];
    iauPom00(xp,yp,sp,as);
    
    // compare the two 3x3 matrices
    int equal = true;
    int different_components = 0;
    double valueAtMaxError;
    double maxd = std::numeric_limits<double>::min();
    for (int r=0; r<3; r++) {
      for (int c=0; c<3; c++) {
        if (!approx_equal(as[r][c], am(r,c))) {
          const double d = std::abs(as[r][c]-am(r,c));
          if (d > maxd) { 
            maxd = d;
            valueAtMaxError = as[r][c];
          }
          equal = false;
          ++different_components;
        }
      }
    }

    if (!equal) {
      ++fails;
      fprintf(stderr, "%23d %20.6e %+.6e\n", different_components, maxd,
              valueAtMaxError);
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
