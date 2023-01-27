#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"fapa03", "fane03", "faur03", "fasa03", "faju03",
                       "fama03", "fae03",  "fave03", "fame03", "faom03",
                       "fad03",  "faf03",  "falp03", "fal03",  "era00", "ee06a",
                       "gmst06", "gmst00"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

bool approx_equal_or_fail(double a, double b, const char *func, int i,
                          double max_allowed = 1e-16) {
  if (!approx_equal(a, b)) {
    fprintf(stderr, "Failing for function: %s at test %d/%d\n", func, i,
            NUM_TESTS);
    fprintf(stderr, "Arguments: %.15f vs %.15f diff: %.5e\n", a, b, a - b);
    if (std::abs(a - b) < max_allowed) {
      fprintf(stderr, "Going on, difference is too small to care!\n");
    } else {
      return false;
    }
  }
  return true;
}

int main(int argc, char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }

  printf("Unit tests for functions:\n");
  for (int i = 0; i < num_funs; i++)
    printf("\t%s\n", funcs[i]);

  double am, as;
  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();

    // add a few seconds to TT to get a random UT1 date
    const auto ut = add_random_seconds(tt, -60e0, 60e0);

    am = fapa03(tt);
    as = iauFapa03(jc);
    assert(approx_equal_or_fail(am, as, "fapa03", i));

    am = fane03(tt);
    as = iauFane03(jc);
    assert(approx_equal_or_fail(am, as, "fane03", i));

    am = faur03(tt);
    as = iauFaur03(jc);
    assert(approx_equal_or_fail(am, as, "faur03", i));

    am = fasa03(tt);
    as = iauFasa03(jc);
    assert(approx_equal_or_fail(am, as, "fasa03", i));

    am = faju03(tt);
    as = iauFaju03(jc);
    assert(approx_equal_or_fail(am, as, "faju03" , i));

    am = fama03(tt);
    as = iauFama03(jc);
    assert(approx_equal_or_fail(am, as, "fama03", i));

    am = fae03(tt);
    as = iauFae03(jc);
    assert(approx_equal_or_fail(am, as, "fae03", i));

    am = fave03(tt);
    as = iauFave03(jc);
    assert(approx_equal_or_fail(am, as, "fave03", i));

    am = fame03(tt);
    as = iauFame03(jc);
    assert(approx_equal_or_fail(am, as, "fame03", i));

    am = faom03(tt);
    as = iauFaom03(jc);
    assert(approx_equal_or_fail(am, as, "faom03", i));

    am = fad03(tt);
    as = iauFad03(jc);
    assert(approx_equal_or_fail(am, as, "fad03", i));

    am = faf03(tt);
    as = iauFaf03(jc);
    assert(approx_equal_or_fail(am, as, "faf03", i));

    am = falp03(tt);
    as = iauFalp03(jc);
    assert(approx_equal_or_fail(am, as, "falp03", i));

    am = fal03(tt);
    as = iauFal03(jc);
    assert(approx_equal_or_fail(am, as, "fal03", i));

    am = era00(tt);
    auto foo = tt.jd_split<dso::TwoPartDate::JdSplitMethod::DT>();
    as = iauEra00(foo._big, foo._small);
    assert(approx_equal_or_fail(am, as, "era00", i));

    auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();
    auto jdut = ut.jd_split<dso::TwoPartDate::JdSplitMethod::DT>();

    am = gmst00(ut, tt);
    as = iauGmst00(jdut._big, jdut._small, jdtt._big, jdtt._small);
    assert(approx_equal_or_fail(am, as, "gmst00", i));
    
    am = gmst06(ut, tt);
    as = iauGmst06(jdut._big, jdut._small, jdtt._big, jdtt._small);
    assert(approx_equal_or_fail(am, as, "gmst06", i));
    
    //am = ee06a(tt);
    //as = iauEe06a(jdtt._big, jdtt._small);
    //assert(approx_equal_or_fail(am, as, "ee06a", i));
  }

  printf("Program %s : All tests passed!\n", argv[0]);
}
