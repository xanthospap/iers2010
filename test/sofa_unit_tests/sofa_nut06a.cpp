#include "datetime/dtfund.hpp"
#include "iau.hpp"
#include "sofa.h"
#include "sofam.h"
#include <cmath>
#include <limits>
#include <random>

//
// test accuracy of (this) nut06a against the original, SOFA implementation
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

  double max_diff[2] = {mind, mind}, ave_diff[2] = {0e0, 0e0};

  for (int t = 0; t < num_tests; t++) {
    double d1[2], d2[2];

    // random MJD (TT)
    const double mjd_tt = uni_mjd(re);

    // call my implementation
    iers2010::sofa::nut06a(dso::mjd0_jd, mjd_tt, d1[0], d1[1]);

    // call SOFA
    iauNut06a(dso::mjd0_jd, mjd_tt, d2+0, d2+1);

    // max_diff and average
    double df;
    for (int i = 0; i < 2; i++) {
      if ((df = std::abs(d1[i] - d2[i])) > max_diff[i])
        max_diff[i] = df;
      ave_diff[i] += d1[i] - d2[i];
    }
  }

  printf("Checking SOFA and *this implementation of NUT06A function. Diffs:\n");
  const char *names[] = {"dpsi", "deps"};
  for (int i = 0; i < 2; i++) {
    printf("\tDelta %s Max:%.15e Average:%+.15e [radians]\n", names[i],
           max_diff[i], ave_diff[i] / num_tests);
    printf("\tDelta %s Max:%.15e Average:%+.15e [arcseconds]\n", names[i],
           max_diff[i] * DR2AS, (ave_diff[i] / num_tests) * DR2AS);
  }

  return 0;
}
