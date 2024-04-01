#include "iau.hpp"
#include "sofa.h"
#include "geodesy/units.hpp"
#include <cstdio>
#include <cassert>
#include <random>

constexpr const int num_tests = 100'000;
constexpr const double MAX_ARCSEC = 1e-12;

int main() {
  double ssofa;
  double smine;
  std::uniform_real_distribution<double> unif(-dso::DPI/5, dso::DPI/5);
  std::default_random_engine re;

  for (int i=0; i<num_tests; i++) {
    /* random X,Y CIP */
    const double x = unif(re);
    const double y = unif(re);
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random(40587, 66154);
    const double dj1 = mjd.imjd() + dso::MJD0_JD;
    const double dj2 = mjd.fractional_days();
    ssofa = iauS06(dj1, dj2, x, y);
    smine = dso::s06(mjd, x, y);
    assert(dso::rad2sec(std::abs(ssofa-smine)) < MAX_ARCSEC);
  }

  printf("%20s %.3e %.5s\n", "s06", MAX_ARCSEC, "OK");
  return 0;
}
