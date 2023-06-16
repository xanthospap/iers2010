#include "iau.hpp"
#include <cstdlib>

double iers2010::sofa::era00(const dso::TwoPartDate &ut1) noexcept {
  /* days since fundamental epoch */
  const dso::TwoPartDate td(ut1.big() - dso::j2000_mjd, ut1.small());
  const double t = td.small() + td.big();

  /* fractional part of day, minus half a day (we need Julian days in the
   * formula but this fractional part is based on MJD)
   */
  const double f = ut1.small();

  /* earth rotation angle at given UT1
   * we are adding here half a day because the fractional part of Julian day
   * is f + .5
   */
  const double theta =
      dso::anp(iers2010::D2PI *
               (f + (0.7790572732640e0 + .5e0) + 0.00273781191135448e0 * t));

  return theta;
}
