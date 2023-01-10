#include "hardisp.hpp"
#include <datetime/dtcalendar.hpp>
#include "geodesy/units.hpp"

//
// Unit test for TDFRPH
//

int main() {
  dso::datetime<dso::nanoseconds> d1(dso::year(2009), dso::month(6),
                                     dso::day_of_month(25),
                                     dso::nanoseconds(0));
  dso::datetime<dso::nanoseconds> d2(
      dso::year(2009), dso::month(6), dso::day_of_month(25), dso::hours(12),
      dso::minutes(1), dso::nanoseconds(45 * 1'000'000'000L));

  double phase, freq;

  iers2010::Hardisp h;
  h(dso::TwoPartDate(d1), dso::TwoPartDate(d1));
  h.tdfrph(dso::DoodsonNumber({2, 0, 0, 0, 0, 0}), phase, freq);

  printf("Frequency: %.15f Phase: %.15f\n", freq/*dso::days_in_julian_cent*/, dso::rad2deg(phase));
  
  h(dso::TwoPartDate(d2), dso::TwoPartDate(d2));
  h.tdfrph(dso::DoodsonNumber({2, 0, 0, 0, 0, 0}), phase, freq);

  printf("Frequency: %.15f Phase: %.15f\n", freq/*dso::days_in_julian_cent*/, dso::rad2deg(phase));

  return 0;
}
