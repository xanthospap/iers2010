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
constexpr const double MAX_ARCSEC_XY = 1e-9;
constexpr const double MAX_ARCSEC_S = 1e-12;

int main() {
  double xsofa, ysofa;
  double xmine, ymine;
  double fargs[14];

  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random(dso::modified_julian_day(47892),
                                           dso::modified_julian_day(66154));
    const double dj1 = mjd.imjd() + dso::MJD0_JD;
    const double dj2 = mjd.fractional_days().days();

    /* SOFA (X,Y) CIP */
    iauXy06(dj1, dj2, &xsofa, &ysofa);

    /* compute (X,Y) CIP and fundamental arguments */
    dso::xycip06a(mjd, xmine, ymine, fargs);
    {
      const double dx = dso::rad2sec(std::abs(xsofa - xmine));
      assert(dx < MAX_ARCSEC_XY);
      const double dy = dso::rad2sec(std::abs(ysofa - ymine));
      assert(dy < MAX_ARCSEC_XY);
    }

    /* use fundamental arguments to compute s */
    const double ssofa = iauS06(dj1, dj2, xsofa, ysofa);
    const double smine = dso::s06(mjd, xmine, ymine, fargs);
    assert(dso::rad2sec(std::abs(ssofa - smine)) < MAX_ARCSEC_S);
  }

  return 0;
}
