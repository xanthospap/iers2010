#include "iau.hpp"
#include "sofa.h"
#include "geodesy/units.hpp"
#include <cstdio>
#include <cassert>

constexpr const int num_tests = 100'000;
constexpr const double MAX_ARCSEC = 1e-9;

int main() {
  double xsofa, ysofa;
  double xmine,ymine;
  double xminex,yminex;
  for (int i=0; i<num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random(40587, 66154);
    const double dj1 = mjd.imjd() + dso::MJD0_JD;
    const double dj2 = mjd.fractional_days();
    iauXy06(dj1, dj2, &xsofa, &ysofa);
    dso::xycip06a(mjd, xmine, ymine);
    assert(std::abs(xsofa-xmine) < MAX_ARCSEC);
    assert(std::abs(ysofa-ymine) < MAX_ARCSEC);
    dso::extra::xycip06a(mjd, xminex, yminex);
    assert(std::abs(xsofa-xminex) < MAX_ARCSEC);
    assert(std::abs(ysofa-yminex) < MAX_ARCSEC);
  }

  return 0;
}
