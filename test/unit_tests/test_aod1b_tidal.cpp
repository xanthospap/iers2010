#include "aod1b.hpp"
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [AOD1B TIDAL PRODUCT]\n", argv[0]);
    return 1;
  }

  dso::StokesCoeffs cs1(180,180,0e0,0e0), cs2(180,180,0e0,0e0);
  // dso::Aod1bIn aptr;

  if (dso::Aod1bIn::get_tidal_wave_coeffs(argv[1], cs1, cs2, nullptr)) {
    fprintf(stderr, "Failed to parse AOD1B file!\n");
    return 1;
  }

  return 0;
}
