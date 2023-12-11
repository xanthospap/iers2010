#include "iau.hpp"
#include "geodesy/units.hpp"

double dso::gmst(const dso::MjdEpoch &tt, const dso::MjdEpoch &ut1) noexcept {
  /* TT Julian centuries since J2000.0. */
  const double t = tt.jcenturies_sinceJ2000();

  /* Greenwich mean sidereal time, IAU 2006. */
  const double gmst = dso::anp(
      dso::era00(ut1) +
      dso::sec2rad(
          0.014506e0 +
          (4612.156534e0 +
           (1.3915817e0 +
            (-0.00000044e0 + (-0.000029956e0 + (-0.0000000368e0) * t) * t) *
                t) *
               t) *
              t));

  return gmst;
}
