#include "iers2010.hpp"

int iers2010::utils::eop_fundarg(double fmjd, double *fargs) noexcept {
  // Convert the input epoch to Julian centuries of TDB since J2000
  double it;
  double ft = std::modf(fmjd, &it);
  const double t = (it - dso::j2000_mjd) / dso::days_in_julian_cent +
                   ft / dso::days_in_julian_cent;

  // Evaluate the vector of the fundamental arguments
  // farg = [ GMST+pi, el, elp, f, d, om ] at t = fmjd

  // 1. Compute GMST ('')
  const double gmst =
      std::fmod(67310.54841e0 + t * ((8640184.812866e0 + 3155760000e0) +
                                     t * (0.093104e0 + t * (-0.0000062e0))),
                86400e0);

  // 24hours are 24*60*60 seconds, hence gmst is ...
  // gmst * 2π / 86400 [radians]

  // 2. Fundamental arguments -> GMST+pi, l, lp, f, d, om
  fargs[0] = std::fmod(gmst * iers2010::D2PI / 86400e0 + iers2010::DPI,
                       iers2010::D2PI);
  iers2010::fundarg(t, fargs + 1);

  return 0;
}
