#include "geodesy/units.hpp"
#include "iau.hpp"

double dso::gmst82(const dso::MjdEpoch &ut1) noexcept {

  constexpr const double half_day_sec = 86400e0 / 2;

  /* Coefficients of IAU 1982 GMST-UT1 model */
  const double A = 24110.54841e0 -  half_day_sec;
  const double B = 8640184.812866e0;
  const double C = 0.093104e0;
  const double D = -6.2e-6;

  /* Julian centuries since fundamental epoch. */
  const double t = ut1.jcenturies_sinceJ2000();

  /* Fractional part of MJD(UT1), in seconds. */
  const double f = ut1.sec_of_day<dso::seconds>();

  /* GMST at this UT1. */
  return dso::anp(
      dso::hsec2rad(((A + (B + (C + D * t) * t) * t) + f) + half_day_sec));
}
