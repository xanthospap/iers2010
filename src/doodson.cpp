#include <cstring>
#include "datetime/dtcalendar.hpp"
#include "geodesy/geoconst.hpp"
#include "doodson.hpp"
#include "iau.hpp"
#include "fundarg.hpp"

char *dso::DoodsonNumber::str(char *buf,
                              bool use_5s_convention) const noexcept {
  int offset = use_5s_convention ? 5 : 0;
  sprintf(buf, "%1d%2d%2d.%2d%2d%2d", iar[0], iar[1] + offset, iar[2] + offset,
          iar[3] + offset, iar[4] + offset, iar[5] + offset);
  return buf;
}

/// @brief Greenwich Mean Sideral Time in [rad], range [0-2π)
/// @warning Note that this is not consistent with IAU2006[A] and IERS2010, 
/// and does not take into account ERA angle and TT time
/// TODO:
/// @see https://github.com/groops-devs/groops/blob/main/source/base/planets.cpp
/// This is not the angle iers2010::gmst??, it does not take into account
/// TT time! However, it seems that this angle is used to compute Doodson
/// arguments
/// See also the 00README_simulation.txt in COST-G benchmark
double dso::gmst_utc(const dso::TwoPartDate &utc) noexcept {
  // julian centuries, at start of day
  const double tu0 = (utc.big() - dso::j2000_mjd) / dso::days_in_julian_cent;
  const double gmst0 =
      (6e0 / 24 + 41e0 / (24 * 60) + 50.54841e0 / (24 * 60 * 60)) +
      (8640184.812866e0 / (24 * 60 * 60)) * tu0 +
      (0.093104e0 / (24 * 60 * 60)) * tu0 * tu0 +
      (-6.2e-6 / (24 * 60 * 60)) * tu0 * tu0 * tu0;
  const double r = 1.002737909350795e0 + 5.9006e-11*tu0 - 5.9e-15*tu0*tu0;
  return dso::anp(dso::D2PI*(gmst0 + r * utc.small()));
}

/*
double dso::DoodsonNumber::phase(
    const dso::TwoPartDate &tt_mjd,
    const dso::TwoPartDate &ut1_mjd) const noexcept {
  const dso::TwoPartDate tt_jd = tt_mjd.jd_sofa();
  const dso::TwoPartDate ut1_jd = ut1_mjd.jd_sofa();
  // Greenwich mean sidereal time (radians)
  const double gmst = iers2010::sofa::gmst06(ut1_jd._big, ut1_jd._small,
                                             tt_jd._big, tt_jd._small);
  // Fundamental arguments
  double fargs[5], beta[6];
  iers2010::fundarg(tt_mjd.jcenturies_sinceJ2000(), fargs);
  // Doodson arguments (store in beta[0:6)
  dso::fundarg2doodson(fargs, gmst, beta);
  return phase(beta);
}*/
