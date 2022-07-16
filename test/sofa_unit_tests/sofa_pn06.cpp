#include "datetime/dtfund.hpp"
#include "iau.hpp"
#include "sofa.h"
#include "sofam.h"
#include <cmath>
#include <limits>
#include <random>

//
// test accuracy of (this) pn06 against the original, SOFA implementation
//
// WARNING!!
// linking against libsofa_c in neccesary to compile this file
//
// this program does not assert anything; it will only report results!
// 
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

// for dpsi, deps (large values)
constexpr const double min_radl = -2e0 * M_PI;
constexpr const double max_radl = 2e0 * M_PI;

// for dpsi, deps (small values)
constexpr const double min_rads = -M_PI / 1e2;
constexpr const double max_rads = M_PI / 1e2;

std::uniform_real_distribution<double> uni_rad_lg(min_radl, max_radl);
std::uniform_real_distribution<double> uni_rad_sm(min_rads, max_rads);

// for MJD
constexpr const double min_mjd = 44239e0; // 01/01/1980
constexpr const double max_mjd = 59945e0; // 01/01/2023
std::uniform_real_distribution<double> uni_mjd(min_mjd, max_mjd);

std::random_device dev;
std::mt19937 rng(dev());
std::uniform_int_distribution<std::mt19937::result_type> rints(1, 6);

int main() {
  // output matrices names:
  // const char *names[] = {"rb", "rp", "rbp", "rn", "rbpn"};

  // Seed with a real random value, if available
  std::random_device r;
  std::default_random_engine re(r());
  constexpr const double mind = std::numeric_limits<double>::min();

  dso::Vector3 unit({1e0,1e0,1e0});

  dso::Vector3 diffs({mind,mind,mind});
  dso::Vector3 avediffs({0e0,0e0,0e0});

  dso::Vector3 diffs_rb=diffs, diffs_rp=diffs, diffs_rbp=diffs, diffs_rn=diffs, diffs_rbpn=diffs;
  dso::Vector3 avediffs_rb=avediffs, avediffs_rp=avediffs, avediffs_rbp=avediffs, avediffs_rn=avediffs, avediffs_rbpn=avediffs;

  for (int t=0; t<num_tests; t++) {
    // random dpsi, deps
    const double dpsi = (t % rints(rng)) ? uni_rad_lg(re)
                              : uni_rad_sm(re);
    const double deps = (t % rints(rng)) ? uni_rad_lg(re)
                              : uni_rad_sm(re);

    // random MJD (TT)
    const double mjd_tt = uni_mjd(re);

    // output
    double epsa1, epsa2;
    dso::Mat3x3 rb1, rp1, rbp1, rn1, rbpn1;

    // call my implementation
    iers2010::sofa::pn06(dso::mjd0_jd, mjd_tt, dpsi, deps, epsa1, rb1, rp1, rbp1, rn1, rbpn1);

    // call SOFA
    double rb2[3][3], rp2[3][3], rbp2[3][3], rn2[3][3], rbpn2[3][3];
    iauPn06(dso::mjd0_jd, mjd_tt, dpsi, deps, &epsa2, rb2, rp2, rbp2, rn2, rbpn2);

    // max_diff and average PER MATRIX
    double v2[3];
    dso::Vector3 v1;
    // -- RB --
    v1 = rb1 * unit;
    iauRxp(rb2, unit.data, v2);
    for (int i = 0; i < 3; i++) {
      double d;
      if ((d = std::abs(v1(i) - v2[i])) > diffs_rb(i)) {
        diffs_rb(i) = d;
      }
      avediffs_rb(i) += v1(i) - v2[i];
    }
    // -- RP --
    v1 = rp1 * unit;
    iauRxp(rp2, unit.data, v2);
    for (int i = 0; i < 3; i++) {
      double d;
      if ((d = std::abs(v1(i) - v2[i])) > diffs_rp(i)) {
        diffs_rp(i) = d;
      }
      avediffs_rp(i) += v1(i) - v2[i];
    }
    // -- RBP --
    v1 = rbp1 * unit;
    iauRxp(rbp2, unit.data, v2);
    for (int i = 0; i < 3; i++) {
      double d;
      if ((d = std::abs(v1(i) - v2[i])) > diffs_rbp(i)) {
        diffs_rbp(i) = d;
      }
      avediffs_rbp(i) += v1(i) - v2[i];
    }
    // -- RN --
    v1 = rn1 * unit;
    iauRxp(rn2, unit.data, v2);
    for (int i = 0; i < 3; i++) {
      double d;
      if ((d = std::abs(v1(i) - v2[i])) > diffs_rn(i)) {
        diffs_rn(i) = d;
      }
      avediffs_rn(i) += v1(i) - v2[i];
    }
    // -- RBPN --
    v1 = rbpn1 * unit;
    iauRxp(rbpn2, unit.data, v2);
    for (int i = 0; i < 3; i++) {
      double d;
      if ((d = std::abs(v1(i) - v2[i])) > diffs_rbpn(i)) {
        diffs_rbpn(i) = d;
      }
      avediffs_rbpn(i) += v1(i) - v2[i];
    }


  }

  printf("Checking SOFA and *this implementation of PN06  function. Diffs:\n");

  printf("Transforming with matrix RB\n");
  printf("\tDelta X Max:%.15e Average:%+.15e [m]\n", diffs_rb(0),
         avediffs_rb(0) / num_tests);
  printf("\tDelta Y Max:%.15e Average:%+.15e [m]\n", diffs_rb(1),
         avediffs_rb(1) / num_tests);
  printf("\tDelta Z Max:%.15e Average:%+.15e [m]\n", diffs_rb(2),
         avediffs_rb(2) / num_tests);

  printf("Transforming with matrix RP\n");
  printf("\tDelta X Max:%.15e Average:%+.15e [m]\n", diffs_rp(0),
         avediffs_rp(0) / num_tests);
  printf("\tDelta Y Max:%.15e Average:%+.15e [m]\n", diffs_rp(1),
         avediffs_rp(1) / num_tests);
  printf("\tDelta Z Max:%.15e Average:%+.15e [m]\n", diffs_rp(2),
         avediffs_rp(2) / num_tests);
  
  printf("Transforming with matrix RBP\n");
  printf("\tDelta X Max:%.15e Average:%+.15e [m]\n", diffs_rbp(0),
         avediffs_rbp(0) / num_tests);
  printf("\tDelta Y Max:%.15e Average:%+.15e [m]\n", diffs_rbp(1),
         avediffs_rbp(1) / num_tests);
  printf("\tDelta Z Max:%.15e Average:%+.15e [m]\n", diffs_rbp(2),
         avediffs_rbp(2) / num_tests);
  
  printf("Transforming with matrix RN\n");
  printf("\tDelta X Max:%.15e Average:%+.15e [m]\n", diffs_rn(0),
         avediffs_rn(0) / num_tests);
  printf("\tDelta Y Max:%.15e Average:%+.15e [m]\n", diffs_rn(1),
         avediffs_rn(1) / num_tests);
  printf("\tDelta Z Max:%.15e Average:%+.15e [m]\n", diffs_rn(2),
         avediffs_rn(2) / num_tests);
  
  printf("Transforming with matrix RBPN\n");
  printf("\tDelta X Max:%.15e Average:%+.15e [m]\n", diffs_rbpn(0),
         avediffs_rbpn(0) / num_tests);
  printf("\tDelta Y Max:%.15e Average:%+.15e [m]\n", diffs_rbpn(1),
         avediffs_rbpn(1) / num_tests);
  printf("\tDelta Z Max:%.15e Average:%+.15e [m]\n", diffs_rbpn(2),
         avediffs_rbpn(2) / num_tests);
  return 0;
}
