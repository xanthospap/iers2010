#include "iau.hpp"
#include "sofa.h"
#include "geodesy/units.hpp"
#include <cstdio>
#include <cassert>

constexpr const int num_tests = 1'000;
//constexpr const double MAX_ARCSEC = 1e-9;

int main() {
  double xsofa, ysofa;
  double xmine,ymine;
  for (int i=0; i<num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random(40587, 66154);
    const double dj1 = mjd.imjd() + dso::MJD0_JD;
    const double dj2 = mjd.fractional_days();
    iauXy06(dj1, dj2, &xsofa, &ysofa);
    dso::xycip06a(mjd, xmine, ymine);
    // assert(std::abs(sofa-mine) < MAX_ARCSEC);
    printf("dX=%+.9f dY=%+.9f [arcsec]\n", dso::rad2sec(xsofa-xmine), dso::rad2sec(ysofa-ymine));
    printf(" X=%+.9f  Y=%+.9f [arcsec]\n", dso::rad2sec(xsofa), dso::rad2sec(ysofa));
  }

  return 0;
}
