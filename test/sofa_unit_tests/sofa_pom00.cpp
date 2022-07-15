#include "datetime/dtfund.hpp"
#include "iau.hpp"
#include "sofa.h"
#include "sofam.h"
#include <cmath>
#include <limits>
#include <random>

//
// test accuracy of (this) pom00 against the original, SOFA implementation
//
// WARNING!!
// linking against libsofa_c in neccesary to compile this file
//
// this program does not assert anything; it will only report results!
// the program must compare two (3x3) matrices, the results that is from the
// call to pom00. The way we check the results of the calls, is using the two 
// matrices (one for each implementation), to transform the vector [1,1,1].
// Discrepancies are check at the results of these, transformed vectors.
// Aka, 
// mat1 = my_pom(...)
// mat2 = sofa_pom(...)
// Vector v(1e0,1e0,1e0)
// v1 = mat1 * v;
// v2 = mat2 * v;
// Diffs = v1 - v2
// 

constexpr const int num_tests = 5000;

// for xp, yp, sp (large values)
constexpr const double min_radl = -2e0 * M_PI;
constexpr const double max_radl = 2e0 * M_PI;

// for xp, yp, sp (small values)
constexpr const double min_rads = -M_PI / 1e2;
constexpr const double max_rads = M_PI / 1e2;

std::uniform_real_distribution<double> uni_rad_lg(min_radl, max_radl);
std::uniform_real_distribution<double> uni_rad_sm(min_rads, max_rads);

std::random_device dev;
std::mt19937 rng(dev());
std::uniform_int_distribution<std::mt19937::result_type> rints(1, 6);

int main() {
  // Seed with a real random value, if available
  std::random_device r;
  std::default_random_engine re(r());
  constexpr const double mind = std::numeric_limits<double>::min();

  dso::Vector3 unit({1e0,1e0,1e0});
  dso::Vector3 diffs({mind,mind,mind});
  dso::Vector3 avediffs({0e0,0e0,0e0});

  for (int t=0; t<num_tests; t++) {
    // random xp, yp, and s'
    const double xp = (t % rints(rng)) ? uni_rad_lg(re)
                              : uni_rad_sm(re);
    const double yp = (t % rints(rng)) ? uni_rad_lg(re)
                              : uni_rad_sm(re);
    const double s  = (t % rints(rng)) ? uni_rad_lg(re)
                              : uni_rad_sm(re);

    // call my implementation
    const auto mat = iers2010::sofa::pom00(xp,yp,s);
    const auto v1 = mat * unit;

    // call SOFA
    double rpom[3][3], v2[3];
    iauPom00(xp, yp, s, rpom);
    iauRxp(rpom, unit.data, v2);

    // max_diff and average
    for (int i = 0; i < 3; i++) {
      double d;
      if ((d = std::abs(v1(i) - v2[i])) > diffs(i)) {
        diffs(i) = d;
      }
      avediffs(i) += v1(i) - v2[i];
    }
  }

  printf("Checking SOFA and *this implementation of POM00  function. Diffs:\n");
  printf("\tDelta X Max:%.15e Average:%+.15e [m]\n", diffs(0),
         avediffs(0) / num_tests);
  printf("\tDelta Y Max:%.15e Average:%+.15e [m]\n", diffs(1),
         avediffs(1) / num_tests);
  printf("\tDelta Z Max:%.15e Average:%+.15e [m]\n", diffs(2),
         avediffs(2) / num_tests);

  return 0;
}
