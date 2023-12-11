#include "eop.hpp"
#include <cstdio>
#include <datetime/calendar.hpp>
#include <datetime/dtfund.hpp>

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
  
  
  dso::MjdEpoch t(57752);
  while (t < dso::MjdEpoch(57754)) {
    dso::EopRecord ie;
    eop.interpolate(t, ie);
    printf("%.12f %.12e %.12e %.12e %.12e %.12e %.12e\n", t.as_mjd(), ie.xp(), ie.yp(), ie.dut(), ie.lod(), ie.dX(), ie.dY());
    t.add_seconds(30e0);
  }

  for (const auto &e : eop.give_me_the_vector()) {
    if (e.t() >= dso::MjdEpoch(57752) && e.t() <= dso::MjdEpoch(57754, dso::FractionalSeconds{3600})) {
      printf("R %.12f %.12e %.12e %.12e %.12e %.12e %.12e\n", e.t().as_mjd(),
             e.xp(), e.yp(), e.dut(), e.lod(), e.dX(), e.dY());
    }   
  }

  return 0;
}
