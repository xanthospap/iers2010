#include "pole_tide.hpp"
#include <cstdio>
#include <cassert>

int main(int argc, char *argv[]) {
  if (argc !=1) {
    fprintf(stderr, "Usage: %s [desaiscopolecoef.txt]\n", argv[0]);
    return 1;
  }

  dso::oceanPoleTide ocp_tide;
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> A_real,
      A_imag, B_real, B_imag;

  ocp_tide::parse_desai02_coeffs(argv[1], A_real, A_imag, B_real, B_imag);

  assert(A_real(0,0)==1e0 
