#include "icgemio.hpp"
#include <cstdio>

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [ICGEM file]\n", argv[0]);
    return 1;
  }

  Icgem gfc(argv[1]);
  StokesCoeffs stokes(80, 60, 0e0, 0e0);
  const datetime<nanoseconds> t(year(2023), month(1), day_of_month(1),
                                nanoseconds(0));

  if (gfc.parse_data(80, 60, t, stokes)) {
    fprintf(stderr, "TEST failed!\n");
    return 1;
  }

  return 0;
}
