#include "iau.hpp"

void iers2010::sofa::pfw06(const dso::TwoPartDate mjd_tt, double &gamb,
                           double &phib, double &psib, double &epsa) noexcept {
  const double t = mjd_tt.jcenturies_sinceJ2000();

  /* P03 bias+precession angles. */
  gamb = dso::sec2rad(
      -0.052928e0 +
      (10.556378e0 +
       (0.4932044e0 +
        (-0.00031238e0 + (-0.000002788e0 + (0.0000000260e0) * t) * t) * t) *
           t) *
          t);

  phib = dso::sec2rad(
      84381.412819e0 +
      (-46.811016e0 +
       (0.0511268e0 +
        (0.00053289e0 + (-0.000000440e0 + (-0.0000000176e0) * t) * t) * t) *
           t) *
          t);

  psib = dso::sec2rad(
      -0.041775e0 +
      (5038.481484e0 +
       (1.5584175e0 +
        (-0.00018522e0 + (-0.000026452e0 + (-0.0000000148e0) * t) * t) * t) *
           t) *
          t);

  epsa = iers2010::sofa::obl06(mjd_tt);

  return;
}
