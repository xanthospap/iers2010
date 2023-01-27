#include "iau.hpp"

void iers2010::sofa::pfw06(const dso::TwoPartDate mjd_tt, double &gamb, double &phib,
                  double &psib, double &epsa) noexcept{
  // double t = ((date1 - dso::j2000_jd) + date2) / dso::days_in_julian_cent;
  const double t = mjd_tt.jcenturies_sinceJ2000();

  // P03 bias+precession angles.
  gamb = (-0.052928 +
          (10.556378 +
           (0.4932044 +
            (-0.00031238 + (-0.000002788 + (0.0000000260) * t) * t) * t) *
               t) *
              t) *
         iers2010::DAS2R;

  phib = (84381.412819 +
          (-46.811016 +
           (0.0511268 +
            (0.00053289 + (-0.000000440 + (-0.0000000176) * t) * t) * t) *
               t) *
              t) *
         iers2010::DAS2R;

  psib = (-0.041775 +
          (5038.481484 +
           (1.5584175 +
            (-0.00018522 + (-0.000026452 + (-0.0000000148) * t) * t) * t) *
               t) *
              t) *
         iers2010::DAS2R;

  epsa = iers2010::sofa::obl06(mjd_tt);

  return;
}
