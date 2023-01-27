#include "iau.hpp"

double iers2010::sofa::era00(const dso::TwoPartDate &mjd_ut1/*double dj1, double dj2*/) noexcept {

  // days since fundamental epoch
  const double t = mjd_ut1._small + (mjd_ut1._big - dso::j2000_jd);

  // fractional part of T in days
  const double f =
      std::fmod(mjd_ut1._big, 1e0) + std::fmod(mjd_ut1._small, 1e0);

  // earth rotation angle at given UT1
  const double theta = dso::anp(
      iers2010::D2PI * (f + 0.7790572732640e0 + 0.00273781191135448e0 * t));

  return theta;
}
