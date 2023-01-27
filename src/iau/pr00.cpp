#include "iau.hpp"

void iers2010::sofa::pr00(const dso::TwoPartDate &mjd_tt, double &dpsipr,
                          double &depspr) noexcept {
  // Precession and obliquity corrections (radians per century)
  constexpr double PRECOR = -0.29965e0 * iers2010::DAS2R;
  constexpr double OBLCOR = -0.02524e0 * iers2010::DAS2R;

  // Interval between fundamental epoch J2000.0 and given date (JC)
  // const double t = ((date1 - dso::j2000_jd) + date2) / dso::days_in_julian_cent;
  const double t = mjd_tt.jcenturies_sinceJ2000();

  // Precession rate contributions with respect to IAU 1976/80.
  dpsipr = PRECOR * t;
  depspr = OBLCOR * t;
}
