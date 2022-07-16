#include "datetime/dtfund.hpp"
#include "iau.hpp"
#include "sofa.h"
#include "sofam.h"
#include <cmath>
#include <limits>
#include <random>

//
// test accuracy of (this) era00 against the original, SOFA implementation
//
// WARNING!!
// linking against libsofa_c in neccesary to compile this file
//
// this program does not assert anything; it will only report results!
//

constexpr const int num_tests = 5000;

// for MJD
constexpr const double min_mjd = 44239e0; // 01/01/1980
constexpr const double max_mjd = 59945e0; // 01/01/2023
std::uniform_real_distribution<double> uni_mjd(min_mjd, max_mjd);

int main() {
  // Seed with a real random value, if available
  std::random_device r;
  std::default_random_engine re(r());
  constexpr const double mind = std::numeric_limits<double>::min();

  double max_diff(mind), ave_diff(0e0);

  for (int t=0; t<num_tests; t++) {
    // random MJD (UT1)
    const double mjd_ut1 = uni_mjd(re);

    // call my implementation
    const double angle1 = iers2010::sofa::era00(dso::mjd0_jd, mjd_ut1);

    // call SOFA
    const double angle2 = iauEra00(dso::mjd0_jd, mjd_ut1);

    // max_diff and average
    double df;
    if ((df=std::abs(angle1-angle2)) > max_diff) {
      max_diff = df;
    }
    ave_diff += angle1 - angle2;
  }

  printf("Checking SOFA and *this implementation of ERA00 function. Diffs:\n");
  printf("\tDelta ERA Max:%.15e Average:%+.15e [radians]\n", max_diff, ave_diff / num_tests);
  printf("\tDelta ERA Max:%.15e Average:%+.15e [arcseconds]\n",
         max_diff * DR2AS, (ave_diff / num_tests) * DR2AS);

return 0;
}
