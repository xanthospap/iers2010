#include "iau.hpp"

double dso::era00(const dso::MjdEpoch &ut1) noexcept {
  /* days since fundamental epoch () */
  constexpr const auto j2000_mjd = dso::MjdEpoch::j2000_mjd();
  const double tu =
      ut1.diff<dso::DateTimeDifferenceType::FractionalDays>(j2000_mjd).days();

  /* fractional part of day */
  const double f = ut1.fractional_days().days();

  /* earth rotation angle at given UT1 */
  const double theta = dso::anp(
      dso::D2PI * (f + 0.7790572732640 + 0.00273781191135448 * tu) + dso::DPI);

  return theta;
}
