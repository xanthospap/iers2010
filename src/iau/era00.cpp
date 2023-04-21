#include "iau.hpp"
#include <cstdlib>
#include <datetime/dtfund.hpp>

double iers2010::sofa::era00(const dso::TwoPartDate &mjd_ut1) noexcept {

  // normalize sot that we are sure the big part is the day, and the small is
  // the fraction
  const dso::TwoPartDate ut1 = mjd_ut1.normalized();

  // days since fundamental epoch
  const double t = ut1.small() + (ut1.big() - dso::j2000_mjd);

  // !Julian! day fraction
  const double f = std::fmod(ut1.small() + .5e0, 1e0);

  // earth rotation angle at given UT1
  const double theta = dso::anp(
      iers2010::D2PI * (f + 0.7790572732640e0 + 0.00273781191135448e0 * t));

  return theta;
}

/*
double iers2010::sofa::era00(const dso::TwoPartDate &mjd_ut1) noexcept {

  // normalize sot that we are sure the big part is the day, and the small is
  // the fraction
  const dso::TwoPartDate jd_ut1(mjd_ut1._big + dso::mjd0_jd, mjd_ut1._small);
  const dso::TwoPartDate ut1 = jd_ut1.normalized();

  // days since fundamental epoch
  const double t = ut1._small + (ut1._big - dso::j2000_jd);

  // !Julian! day fraction
  const double f = std::fmod(ut1._big, 1.0) + std::fmod(ut1._small, 1.0);

  // earth rotation angle at given UT1
  const double theta = dso::anp(
      iers2010::D2PI * (f + 0.7790572732640e0 + 0.00273781191135448e0 * t));

  return theta;
}
*/
