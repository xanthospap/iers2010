#include "iau.hpp"
#include "sofa.h"
#include <random>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

constexpr const int num_tests = 100'000;
constexpr const double MAX_ARCSEC = 1e-7;

std::uniform_real_distribution<double> unif(-60e0, 60e0);
std::default_random_engine re;

int main() {
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch (TT) */
    const auto mjd_tt = dso::MjdEpoch::random(dso::modified_julian_day(47892),
                                              dso::modified_julian_day(66154));
    const double dj1_tt = mjd_tt.imjd() + dso::MJD0_JD;
    const double dj2_tt = mjd_tt.fractional_days().days();

    /* UT1 to assign gmst computation +- 1min from TT epoch */
    const auto mjd_ut1 = mjd_tt.add_seconds(dso::FractionalSeconds(unif(re)));
    const double dj1_ut1 = mjd_ut1.imjd() + dso::MJD0_JD;
    const double dj2_ut1 = mjd_ut1.fractional_days().days();

    const double sofa = iauGmst06(dj1_ut1, dj2_ut1, dj1_tt, dj2_tt);
    const double mine = dso::gmst(mjd_tt, mjd_ut1);

    assert(dso::rad2sec(std::abs(sofa - mine)) < MAX_ARCSEC);

    const double sofa_wrong = iauGmst00(dj1_ut1, dj2_ut1, dj1_tt, dj2_tt);
    assert(dso::rad2sec(std::abs(sofa_wrong - mine)) >
           dso::rad2sec(std::abs(sofa - mine)));
  }

  return 0;
}
