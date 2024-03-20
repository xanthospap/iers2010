#include "iau.hpp"
#include "geodesy/units.hpp"

double dso::dgmst(const dso::MjdEpoch &tt) noexcept {
  const double t = tt.jcenturies_sinceJ2000();
  const double b =
      dso::sec2rad(
          4612.156534e0 +
          (1.3915817e0 +
           (-0.00000044e0 + (-0.000029956e0 + (-0.0000000368e0) * t) * t) * t) *
              t) /
      36525e0;
  constexpr const double t2pi = dso::detail::AngleUnitTraits<
      dso::detail::AngleUnit::Radians>::full_circle();
  const double a = t2pi * 2e0 * 1.00273781191135448;
  return a + b;
}
