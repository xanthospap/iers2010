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

int main(int argc, char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }

  printf("***************************************************************\n");
  printf("Unit tests for functions:\n");
  for (int i = 0; i < num_funs; i++)
    printf("\t%s\n", funcs[i]);
  printf("---------------------------------------------------------------\n");

  double am, as;
  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fapa03(tt);
    as = iauFapa03(jc);
    assert(approx_equal(am, as, "fapa03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fane03(tt);
    as = iauFane03(jc);
    assert(approx_equal(am, as, "fane03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = faur03(tt);
    as = iauFaur03(jc);
    assert(approx_equal(am, as, "faur03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fasa03(tt);
    as = iauFasa03(jc);
    assert(approx_equal(am, as, "fasa03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = faju03(tt);
    as = iauFaju03(jc);
    assert(approx_equal(am, as, "faju03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fama03(tt);
    as = iauFama03(jc);
    assert(approx_equal(am, as, "fama03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fae03(tt);
    as = iauFae03(jc);
    assert(approx_equal(am, as, "fae03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fave03(tt);
    as = iauFave03(jc);
    assert(approx_equal(am, as, "fave03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fame03(tt);
    as = iauFame03(jc);
    assert(approx_equal(am, as, "fame03"));
  }
  
  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = faom03(tt);
    as = iauFaom03(jc);
    assert(approx_equal(am, as, "faom03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fad03(tt);
    as = iauFad03(jc);
    assert(approx_equal(am, as, "fad03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = faf03(tt);
    as = iauFaf03(jc);
    assert(approx_equal(am, as, "faf03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = falp03(tt);
    as = iauFalp03(jc);
    assert(approx_equal(am, as, "falp03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fal03(tt);
    as = iauFal03(jc);
    assert(approx_equal(am, as, "fal03"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    am = era00(tt);
    auto foo = tt.jd_split<dso::TwoPartDate::JdSplitMethod::DT>();
    as = iauEra00(foo._big, foo._small);
    assert(approx_equal(am, as, "era00"));
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    // add a few seconds to TT to get a random UT1 date
    const auto ut = add_random_seconds(tt, -60e0, 60e0);
    auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();
    auto jdut = ut.jd_split<dso::TwoPartDate::JdSplitMethod::DT>();

    am = gmst00(ut, tt);
    as = iauGmst00(jdut._big, jdut._small, jdtt._big, jdtt._small);
    assert(approx_equal(am, as, "gmst00"));
  }
    
  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    // add a few seconds to TT to get a random UT1 date
    const auto ut = add_random_seconds(tt, -60e0, 60e0);
    auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();
    auto jdut = ut.jd_split<dso::TwoPartDate::JdSplitMethod::DT>();

    am = gmst06(ut, tt);
    as = iauGmst06(jdut._big, jdut._small, jdtt._big, jdtt._small);
    assert(approx_equal(am, as, "gmst06"));
  }
    
  for (int i = 0; i < NUM_TESTS; i++) {
    // random date (MJD, TT)
    const auto tt = random_mjd();
    // add a few seconds to TT to get a random UT1 date
    auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();

    am = ee06a(tt);
    as = iauEe06a(jdtt._big, jdtt._small);
    assert(approx_equal(am, as, "ee06a"));
  }

  printf("Program %s : All tests passed!\n", argv[0]);
  printf("***************************************************************\n");
}
