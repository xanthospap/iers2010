#include "iau.hpp"
#include "iersc.hpp"
#include <algorithm>
#include <cmath>

double iers2010::sofa::era00(double dj1, double dj2) noexcept {
  double d1 = dj1;
  double d2 = dj2;
  if (dj1 >= dj2)
    std::swap(d1, d2);

  // days since fundamental epoch
  double t = d1 + (d2 - iers2010::DJ00);

  // fractional part of T in days
  double f = std::fmod(d1, 1e0) + std::fmod(d2, 1e0);

  // earth rotation angle at given UT1
  double theta = iers2010::nang_02pi(
      iers2010::D2PI * (f + 0.7790572732640e0 + 0.00273781191135448e0 * t));

  return theta;
}