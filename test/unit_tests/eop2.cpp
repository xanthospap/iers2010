#include "eop.hpp"
#include <cassert>
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage %s [EOP C04 FILE]\n", argv[0]);
    return 1;
  }

  dso::MjdEpoch t1(57750);
  dso::MjdEpoch t2(57759);
  dso::EopSeries eop;

  if (dso::parse_iers_C04(argv[1], t1, t2, eop)) {
    fprintf(stderr, "ERROR Failed parsing eop file\n");
    return 1;
  }
  assert(eop.num_entries() == 9);

  dso::MjdEpoch t(57750);
  dso::EopRecord ie;
  const auto &veop = eop.give_me_the_vector();

  /* interpolate, prior to first date in series */
  assert(eop.interpolate(t, ie) ==
         dso::EopSeries::EopInterpolationResult::OutOfBoundsPrior);
  assert(ie.t() == t);
  assert(ie.xp() == veop.begin()->xp());
  assert(ie.yp() == veop.begin()->yp());
  assert(ie.dut() == veop.begin()->dut());
  assert(ie.lod() == veop.begin()->lod());
  assert(ie.dX() == veop.begin()->dX());
  assert(ie.dY() == veop.begin()->dY());

  /* interpolate later than the last date in series */
  t = dso::MjdEpoch(57759, dso::FractionalSeconds{60});
  assert(eop.interpolate(t, ie) ==
         dso::EopSeries::EopInterpolationResult::OutOfBoundsLater);
  assert(ie.t() == t);
  assert(ie.xp() == (veop.end() - 1)->xp());
  assert(ie.yp() == (veop.end() - 1)->yp());
  assert(ie.dut() == (veop.end() - 1)->dut());
  assert(ie.lod() == (veop.end() - 1)->lod());
  assert(ie.dX() == (veop.end() - 1)->dX());
  assert(ie.dY() == (veop.end() - 1)->dY());

  /* only linear interpolation possible */
  t = dso::MjdEpoch(57750, dso::FractionalSeconds{60});
  assert(eop.interpolate(t, ie) ==
         dso::EopSeries::EopInterpolationResult::Linear);
  assert(ie.t() == t);
  t = dso::MjdEpoch(57758, dso::FractionalSeconds{0});
  assert(eop.interpolate(t, ie) ==
         dso::EopSeries::EopInterpolationResult::Linear);
  assert(ie.t() == t);

  /* maximum degree that can be achieved = 3 */
  t = dso::MjdEpoch(57751, dso::FractionalSeconds{60});
  assert(eop.interpolate(t, ie) ==
         dso::EopSeries::EopInterpolationResult::PolyDegreeDescreased);
  assert(ie.t() == t);
  t = dso::MjdEpoch(57757, dso::FractionalSeconds{0});
  assert(eop.interpolate(t, ie) ==
         dso::EopSeries::EopInterpolationResult::PolyDegreeDescreased);
  assert(ie.t() == t);

  /* maximum degree that can be achieved = 5 */
  t = dso::MjdEpoch(57752, dso::FractionalSeconds{60});
  assert(eop.interpolate(t, ie) ==
         dso::EopSeries::EopInterpolationResult::PolyDegreeRequested);
  assert(ie.t() == t);
  t = dso::MjdEpoch(57756, dso::FractionalSeconds{0});
  assert(eop.interpolate(t, ie) ==
         dso::EopSeries::EopInterpolationResult::PolyDegreeRequested);
  assert(ie.t() == t);
  
  return 0;
}
