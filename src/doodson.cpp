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
