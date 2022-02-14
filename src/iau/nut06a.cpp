#include "iau.hpp"
#include "iersc.hpp"

void iers2010::sofa::nut06a(double date1, double date2, double &dpsi,
                            double &deps) noexcept {
  // Interval between fundamental date J2000.0 and given date (JC).
  const double t = ((date1 - iers2010::DJ00) + date2) / iers2010::DJC;

  // Factor correcting for secular variation of J2.
  const double fj2 = -2.7774e-6 * t;

  // Obtain IAU 2000A nutation.
  double dp, de;
  iers2010::sofa::nut00a(date1, date2, dp, de);

  // Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5).
  dpsi = dp + dp * (0.4697e-6 + fj2);
  deps = de + de * fj2;

  return;
}