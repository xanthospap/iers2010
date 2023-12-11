#include "eop.hpp"
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage %s [EOP C04 FILE]\n", argv[0]);
    return 1;
  }

  dso::MjdEpoch t1(57753);
  dso::MjdEpoch t2(57755);
  dso::EopSeries eop;

  if (dso::parse_iers_C04(argv[1], t1, t2, eop)) {
    fprintf(stderr, "ERROR Failed parsing eop file\n");
    return 1;
  }

  assert(eop.num_entries() == 2);

  return 0;
}
