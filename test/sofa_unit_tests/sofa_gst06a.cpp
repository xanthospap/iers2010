#include "datetime/dtfund.hpp"
#include "iau.hpp"
#include "sofa.h"
#include "sofam.h"
#include <cmath>
#include <limits>
#include <random>

//
// test accuracy of (this) gst06a against the original, SOFA implementation
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

// for TT - UT1 difference
constexpr const double min_sec = -59.9e0;
constexpr const double max_sec = 59.9e0;
std::uniform_real_distribution<double> uni_sec(min_sec, max_sec);

int main() {
  // Seed with a real random value, if available
  std::random_device r;
  std::default_random_engine re(r());
  constexpr const double mind = std::numeric_limits<double>::min();

  double max_diff(mind), ave_diff(0e0);

  for (int t=0; t<num_tests; t++) {
    // random MJD (TT)
    const double mjd_tt = uni_mjd(re);

    // random TT - UT1 difference in seconds (range -59.9 sec to  59.9 sec)
    const double dsec = uni_sec(re);

    // construct a UT1 date
    const double mjd_ut1 = mjd_tt + (dsec / 86400e0);

    // call my implementation
    const double g06a1 = iers2010::sofa::gst06a(dso::mjd0_jd, mjd_ut1, dso::mjd0_jd, mjd_tt);

    // call SOFA
    const double g06a2 = iauGst06a(dso::mjd0_jd, mjd_ut1, dso::mjd0_jd, mjd_tt);

    // max_diff and average
    double df;
    if ((df=std::abs(g06a1-g06a2)) > max_diff) max_diff = df;
    ave_diff += g06a1 - g06a2;
  }

  printf("Checking SOFA and *this implementation of GST06A function. Diffs:\n");
  printf("\tDelta GST Max:%.15e Average:%+.15e [radians]\n", max_diff, ave_diff / num_tests);
  printf("\tDelta GST Max:%.15e Average:%+.15e [arcseconds]\n",
         max_diff * DR2AS, (ave_diff / num_tests) * DR2AS);

return 0;
}
