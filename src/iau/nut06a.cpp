#include "iau.hpp"

void iers2010::sofa::nut06a(const dso::TwoPartDate &mjd_tt, double &dpsi,
                            double &deps) noexcept {
  /* Interval between fundamental date J2000.0 and given date (JC). */
  const double t = mjd_tt.jcenturies_sinceJ2000();

  /* Factor correcting for secular variation of J2. */
  const double fj2 = -2.7774e-6 * t;

  /* Obtain IAU 2000A nutation. */
  double dp, de;
  iers2010::sofa::nut00a(mjd_tt, dp, de);

  /* Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5).*/
  dpsi = dp + dp * (0.4697e-6 + fj2);
  deps = de + de * fj2;

  return;
}
