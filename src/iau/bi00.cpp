#include "iau.hpp"

/* note that this is a const function, all variables are pre-computed */
void iers2010::sofa::bi00(double &dpsibi, double &depsbi,
                          double &dra) noexcept {
  /* The frame bias corrections in longitude and obliquity */
  constexpr const double DPBIAS = dso::sec2rad(-0.041775e0);
  constexpr const double DEBIAS = dso::sec2rad(-0.0068192e0);

  /* The ICRS RA of the J2000.0 equinox (Chapront et al., 2002) */
  constexpr const double DRA0 = dso::sec2rad(-0.0146e0);

  /* Return the results (which are fixed). */
  dpsibi = DPBIAS;
  depsbi = DEBIAS;
  dra = DRA0;
}
