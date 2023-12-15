#include "iau.hpp"
#include "sofa.h"
#include "geodesy/units.hpp"
#include <cstdio>
#include <cassert>
#include <limits>

constexpr const int num_tests = 100'000;
constexpr const double MAX_ARCSEC = 1e-9;
constexpr const double MIND = std::numeric_limits<double>::min();

int main() {
  double xsofa, ysofa;
  double xmine,ymine;
  double xminex,yminex;
  double maxsofav1[2] = {MIND, MIND};
  double maxsofav2[2] = {MIND, MIND};
  double maxv1v2[2] = {MIND, MIND};

  for (int i=0; i<num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random(40587, 66154);
    const double dj1 = mjd.imjd() + dso::MJD0_JD;
    const double dj2 = mjd.fractional_days();
    /* SOFA */
    iauXy06(dj1, dj2, &xsofa, &ysofa);
    /* this implementation: v1 i.e. default */
    dso::xycip06a(mjd, xmine, ymine);
    {
      const double dx = std::abs(xsofa-xmine);
      if (maxsofav1[0] <= dx) maxsofav1[0] = dx;
      assert(dx < MAX_ARCSEC);
      const double dy = std::abs(ysofa-ymine);
      if (maxsofav1[1] <= dy) maxsofav1[1] = dy;
      assert(dy < MAX_ARCSEC);
    }
    /* this implementation: v2 i.e. seperate x&y series */
    dso::extra::xycip06a(mjd, xminex, yminex);
    {
      const double dx = std::abs(xsofa-xminex);
      if (maxsofav2[0] <= dx) maxsofav2[0] = dx;
      assert(dx < MAX_ARCSEC);
      const double dy = std::abs(ysofa-yminex);
      if (maxsofav1[1] <= dy) maxsofav1[1] = dy;
      assert(dy < MAX_ARCSEC);
    }
    /* compare v1 to v2 */
    {
      const double dx = std::abs(xmine-xminex);
      if (maxv1v2[0] <= dx) maxv1v2[0] = dx;
      assert(dx < MAX_ARCSEC);
      const double dy = std::abs(ymine-yminex);
      if (maxv1v2[1] <= dy) maxv1v2[1] = dy;
      assert(dy < MAX_ARCSEC);
    }
  }

  printf("Max differences between implementations in [microarcsec]\n");
  printf("SOFA - v1: %.5e %.5e\n", maxsofav1[0], maxsofav1[1]);
  printf("SOFA - v2: %.5e %.5e\n", maxsofav2[0], maxsofav2[1]);
  printf("v1 - v2  : %.5e %.5e\n", maxv1v2[0], maxv1v2[1]);

  return 0;
}
