#include "iau.hpp"

void iers2010::sofa::bi00(double &dpsibi, double &depsbi,
                          double &dra) noexcept {
  // The frame bias corrections in longitude and obliquity
  constexpr double DPBIAS = -0.041775e0 * DAS2R;
  constexpr double DEBIAS = -0.0068192e0 * DAS2R;

  // The ICRS RA of the J2000.0 equinox (Chapront et al., 2002)
  constexpr double DRA0 = -0.0146e0 * DAS2R;

  // Return the results (which are fixed).
  dpsibi = DPBIAS;
  depsbi = DEBIAS;
  dra = DRA0;
}
