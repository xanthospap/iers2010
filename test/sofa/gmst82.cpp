#include "iau.hpp"
#include "sofa.h"
#include <cstdio>
#include <cassert>

constexpr const int num_tests = 100'000;
constexpr const double MAX_ARCSEC = 1e-7;

int main() {
  for (int i=0; i<num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random(40587, 66154);
    const double dj1 = mjd.imjd() + dso::MJD0_JD;
    const double dj2 = mjd.fractional_days();
    const double sofa = iauGmst82(dj1, dj2);
    const double mine = dso::gmst82(mjd);
    assert(dso::rad2sec(std::abs(sofa-mine)) < MAX_ARCSEC);
  }

  printf("%.20s %.3e %.5s\n", "gmst82", MAX_ARCSEC, "OK");
  return 0;
}
