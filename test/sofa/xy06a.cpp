#include "geodesy/units.hpp"
#include "iau.hpp"
#include "sofa.h"
#include <cassert>
#include <limits>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

constexpr const int num_tests = 10'000;
constexpr const double MAX_ARCSEC_V1 = 1e-9;
constexpr const double MAX_ARCSEC_V2 = 1e-8;

int main()
{
  double xsofa, ysofa;
  double xmine, ymine;
  double xminex, yminex;

  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random(dso::modified_julian_day(47892),
        dso::modified_julian_day(66154));
    const double dj1 = mjd.imjd() + dso::MJD0_JD;
    const double dj2 = mjd.fractional_days().days();

    /* SOFA */
    iauXy06(dj1, dj2, &xsofa, &ysofa);

    /* this implementation: v1 i.e. default */
    dso::xycip06a(mjd, xmine, ymine);
    {
      const double dx = dso::rad2sec(std::abs(xsofa - xmine));
      assert(dx < MAX_ARCSEC_V1);
      const double dy = dso::rad2sec(std::abs(ysofa - ymine));
      assert(dy < MAX_ARCSEC_V1);
    }

    /* this implementation: v2 i.e. seperate x&y series */
    dso::extra::xycip06a(mjd, xminex, yminex);
    {
      const double dx = dso::rad2sec(std::abs(xsofa - xminex));
      assert(dx < MAX_ARCSEC_V2);
      const double dy = dso::rad2sec(std::abs(ysofa - yminex));
      assert(dy < MAX_ARCSEC_V2);
    }

  }

  return 0;
}
