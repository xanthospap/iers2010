#include "icgemio.hpp"
#include <cstdio>

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [ICGEM file]\n", argv[0]);
    return 1;
  }

  Icgem gfc(argv[1]);
  StokesCoeffs stokes(0, 0, 0e0, 0e0);

  int size = 120;
  for (int iy = 1995; iy < 2025; iy++) {
    /* set request date */
    const datetime<nanoseconds> t(year(iy), month(1), day_of_month(1),
                                  nanoseconds(0));
    /* parse data for the give date and degree/order */
    if (gfc.parse_data(size, size - 1, t, stokes)) {
      fprintf(stderr, "TEST failed!\n");
      return 1;
    }

    /* print results for C(3,2) and S(3,2) */
    printf("C = %+.15e S = %+.15e\n", stokes.C(3, 2), stokes.S(3, 2));

    /* check the coefficients size */
    assert(stokes.max_degree() == size);
    assert(stokes.max_order() == size-1);

    /* decrease parsing size */
    size -= 5;
    if (size < 3) size = 120;
  }

  return 0;
}
