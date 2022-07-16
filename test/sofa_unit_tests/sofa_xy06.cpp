#include "datetime/dtfund.hpp"
#include "iau.hpp"
#include "sofa.h"
#include "sofam.h"
#include <cmath>
#include <limits>
#include <random>

//
// test accuracy of (this) xy06 against the original, SOFA implementation
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

  double max_diff[] = {mind,mind};
  double ave_diff[] = {0e0, 0e0};
  double xy1[2], xy2[2];

  for (int t=0; t<num_tests; t++) {
    // random MJD (TT)
    const double mjd_tt = uni_mjd(re);

    // call my implementation
    iers2010::sofa::xy06(dso::mjd0_jd, mjd_tt, xy1[0], xy1[1]);

    // call SOFA
    iauXy06(dso::mjd0_jd, mjd_tt, xy2+0, xy2+1);

    // max_diff and average
    for (int i = 0; i < 2; i++) {
      double d;
      if ((d = std::abs(xy1[i] - xy2[i])) > max_diff[i]) {
        max_diff[i] = d;
      }
      ave_diff[i] += xy1[i] - xy2[i];
    }
  }

  printf("Checking SOFA and *this implementation of XY06 function. Diffs:\n");
  printf("\tDelta X Max:%.15e Average:%+.15e \n", max_diff[0],
         ave_diff[0] / num_tests);
  printf("\tDelta Y Max:%.15e Average:%+.15e \n", max_diff[1],
         ave_diff[1] / num_tests);
  return 0;
}
