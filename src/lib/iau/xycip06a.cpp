#include "iau.hpp"
#include "fundarg.hpp"
#include "geodesy/units.hpp"
#include <cstring>

int dso::extra::xycip06a(const dso::MjdEpoch &tt, double &xcip,
                  double &ycip, double *outargs) noexcept {
  /* Interval between fundamental date J2000.0 and given date. */
  const double t = tt.jcenturies_sinceJ2000();

  /* Fundamental arguments (IERS 2003) */
  const double fa[14] = {
      dso::iers2010::fal03(t),
      dso::iers2010::falp03(t),
      dso::iers2010::faf03(t),
      dso::iers2010::fad03(t),
      dso::iers2010::faom03(t),
      dso::iers2010::fame03(t),
      dso::iers2010::fave03(t),
      dso::iers2010::fae03(t),
      dso::iers2010::fama03(t),
      dso::iers2010::faju03(t),
      dso::iers2010::fasa03(t),
      dso::iers2010::faur03(t),
      dso::iers2010::fane03(t),
      dso::iers2010::fapa03(t),
  };

  /* X-coordinate of CIP [μas] */
  const double xcip_mas = dso::detail::xcip06a(fa, t);

  /* Y-coordinate of CIP [μas] */
  const double ycip_mas = dso::detail::ycip06a(fa, t);

  /* transform to [rad] and return */
  xcip = dso::sec2rad(xcip_mas*1e-6);
  ycip = dso::sec2rad(ycip_mas*1e-6);

  /* if we are given a large enough array, store the fundamental arguments */
  if (outargs) 
    std::memcpy(outargs, fa, sizeof(double) * 14);

  return 0;
}
