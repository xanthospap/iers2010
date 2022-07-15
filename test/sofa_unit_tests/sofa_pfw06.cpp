#include "datetime/dtfund.hpp"
#include "iau.hpp"
#include "sofa.h"
#include "sofam.h"
#include <cmath>
#include <limits>
#include <random>

//
// test accuracy of (this) pfw06 against the original, SOFA implementation
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

  double max_diff[4];
  for (int i = 0; i < 4; i++)
    max_diff[i] = mind;
  double ave_diff[4];
  for (int i = 0; i < 4; i++)
    ave_diff[i] = 0e0;
  const char *names[] = {"gamb", "phib", "psib", "epsa"};

  for (int t = 0; t < num_tests; t++) {
    // random MJD (TT)
    const double mjd_tt = uni_mjd(re);

    double angles1[4], angles2[4];

    // call my implementation
    iers2010::sofa::pfw06(dso::mjd0_jd, mjd_tt, angles1[0], angles1[1],
                          angles1[2], angles1[3]);

    // call SOFA
    iauPfw06(dso::mjd0_jd, mjd_tt, angles2 + 0, angles2 + 1, angles2 + 2,
             angles2 + 3);

    // max_diff and average
    double df;
    for (int i = 0; i < 4; i++) {
      if ((df = std::abs(angles1[i] - angles2[i])) > max_diff[i])
        max_diff[i] = df;
      ave_diff[i] += angles1[i] - angles2[i];
    }
  }

  printf("Checking SOFA and *this implementation of PFW06 function. Diffs:\n");
  for (int i = 0; i < 4; i++) {
    printf("\tDelta %s Max:%.15e Average:%+.15e [radians]\n", names[i],
           max_diff[i], ave_diff[i] / num_tests);
    printf("\tDelta %s Max:%.15e Average:%+.15e [arcseconds]\n", names[i],
           max_diff[i] * DR2AS, (ave_diff[i] / num_tests) * DR2AS);
  }

  return 0;
}
