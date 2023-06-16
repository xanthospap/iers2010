#include "iau.hpp"

void iers2010::sofa::pr00(const dso::TwoPartDate &mjd_tt, double &dpsipr,
                          double &depspr) noexcept {
  /* Precession and obliquity corrections (radians per century) */
  constexpr const double PRECOR = dso::sec2rad(-0.29965e0);
  constexpr const double OBLCOR = dso::sec2rad(-0.02524e0);

  /* Interval between fundamental epoch J2000.0 and given date (JC) */
  const double t = mjd_tt.jcenturies_sinceJ2000();

  /* Precession rate contributions with respect to IAU 1976/80. */
  dpsipr = PRECOR * t;
  depspr = OBLCOR * t;
}
