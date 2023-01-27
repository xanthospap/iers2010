#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 1000;

const double EPSILON = 1e-16;

const char *funcs[] = {"pr00", "numat", "pnm00a", "pn00a", "pn00", "pn06", "p06e", "bp00"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main(int argc, char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }

  printf("Unit tests for functions:\n");
  for (int i=0; i<num_funs; i++) printf("\t%s\n", funcs[i]);

  double am, as;
  double rot[3][3];
  Eigen::Matrix<double, 3, 3> r1, r2, r3, r4, r5;
  double m1[3][3], m2[3][3], m3[3][3], m4[3][3], m5[3][3];
  double dm[20], ds[20];
  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    // const auto jc = tt.jcenturies_sinceJ2000();

    // add a few seconds to TT to get a random UT1 date
    // const auto ut = add_random_seconds(tt, -60e0, 60e0);

    // split TT JD in a "J2000" fashion
    auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();

    pr00(tt, dm[0],dm[1]);
    iauPr00(jdtt._big, jdtt._small, &ds[0],&ds[1]);
    assert(approx_equal(dm[0], ds[0]/*, EPSILON*/));
    assert(approx_equal(dm[1], ds[1]/*, EPSILON*/));

    dm[0] = ds[0] = random_angle(-dso::DPI/5e0, dso::DPI/5);
    dm[1] = ds[1] = random_angle(-dso::DPI/5e0, dso::DPI/5);
    dm[2] = ds[2] = random_angle(-dso::DPI/5e0, dso::DPI/5);
    iauNumat(ds[0],ds[1],ds[2],rot);
    assert(approx_equal(numat(dm[0],dm[1],dm[2]), rot, EPSILON, true));

    const double dpsi = ds[0];
    const double deps = ds[1];
    pn00(tt, dpsi, deps, am, r1, r2, r3, r4, r5);
    iauPn00(jdtt._big, jdtt._small, dpsi, deps, &as, m1, m2, m3, m4, m5);
    assert(approx_equal(r1, m1, EPSILON));
    assert(approx_equal(r2, m2, EPSILON));
    assert(approx_equal(r3, m3, EPSILON));
    assert(approx_equal(r4, m4, EPSILON));
    assert(approx_equal(r5, m5, EPSILON));
    assert(approx_equal(am, as, EPSILON));
    
    iauPnm00a(jdtt._big, jdtt._small, rot);
    assert(approx_equal(pnm00a(tt), rot, EPSILON, true));

    pn00a(tt, dm[0], dm[1], dm[2], r1, r2, r3, r4, r5);
    iauPn00a(jdtt._big, jdtt._small, &ds[0], &ds[1], &ds[2], m1, m2, m3, m4, m5);
    assert(approx_equal(r1, m1, EPSILON));
    assert(approx_equal(r2, m2, EPSILON));
    assert(approx_equal(r3, m3, EPSILON));
    assert(approx_equal(r4, m4, EPSILON));
    assert(approx_equal(r5, m5, EPSILON));
    assert(approx_equal(dm[0], ds[0], EPSILON));
    assert(approx_equal(dm[1], ds[1], EPSILON));
    assert(approx_equal(dm[2], ds[2], EPSILON));


    pn06(tt, dpsi, deps, am, r1, r2, r3, r4, r5);
    iauPn06(jdtt._big, jdtt._small, dpsi, deps, &as, m1, m2, m3, m4, m5);
    assert(approx_equal(r1, m1, EPSILON));
    assert(approx_equal(r2, m2, EPSILON));
    assert(approx_equal(r3, m3, EPSILON));
    assert(approx_equal(r4, m4, EPSILON));
    assert(approx_equal(r5, m5, EPSILON));
    assert(approx_equal(am, as, EPSILON));

    p06e(tt, dm[0], dm[1], dm[2], dm[3], dm[4], dm[5], dm[6], dm[7], dm[8],
         dm[9], dm[10], dm[11], dm[12], dm[13], dm[14], dm[15]);
    iauP06e(jdtt._big, jdtt._small, &ds[0], &ds[1], &ds[2], &ds[3], &ds[4], &ds[5],
            &ds[6], &ds[7], &ds[8], &ds[9], &ds[10], &ds[11], &ds[12], &ds[13], &ds[14],
            &ds[15]);
    for (int j = 0; j < 16; j++)
      assert(approx_equal(dm[j], ds[j], EPSILON));

    bp00(tt, r1, r2, r3);
    iauBp00(jdtt._big, jdtt._small, m1, m2, m3);
    assert(approx_equal(r1, m1, EPSILON));
    assert(approx_equal(r2, m2, EPSILON));
    assert(approx_equal(r3, m3, EPSILON));
  }
  
  printf("Program %s : All tests passed!\n", argv[0]);
  return 0;
}
