#include "pole_tide.hpp"
#include <random>

using namespace dso;

int main(int argc, char *argv[]) {
  
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [desaiscopolecoef.txt]\n", argv[0]);
    return 1;
  }

  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> A_real_ref(
      pole_tide_details::MAX_DEGREE_DESAI_2002 + 1,
      pole_tide_details::MAX_DEGREE_DESAI_2002 + 1);
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> A_imag_ref(
      pole_tide_details::MAX_DEGREE_DESAI_2002 + 1,
      pole_tide_details::MAX_DEGREE_DESAI_2002 + 1);
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> B_real_ref(
      pole_tide_details::MAX_DEGREE_DESAI_2002 + 1,
      pole_tide_details::MAX_DEGREE_DESAI_2002 + 1);
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> B_imag_ref(
      pole_tide_details::MAX_DEGREE_DESAI_2002 + 1,
      pole_tide_details::MAX_DEGREE_DESAI_2002 + 1);
  
  /* if the test test_parse_desai02_coeffs.cpp has passed, this should work well */
  assert(!pole_tide_details::parse_desai02_coeffs(
      argv[1], pole_tide_details::MAX_DEGREE_DESAI_2002,
      pole_tide_details::MAX_DEGREE_DESAI_2002, A_real_ref, A_imag_ref, B_real_ref,
      B_imag_ref));

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distr(1, 360);

  for (int i = 0; i < 10; i++) {
    const int maxDegree = distr(gen);
    std::uniform_int_distribution<> distr2(1, maxDegree);
    const int maxOrder = distr2(gen);

    /* let the parse_desai02_coeffs do the allocation */
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> A_real(1, 1);
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> A_imag(1, 1);
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> B_real(1, 1);
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> B_imag(1, 1);

    assert(!pole_tide_details::parse_desai02_coeffs(
        argv[1], maxDegree, maxOrder, A_real, A_imag, B_real, B_imag));

    /* validate all */
    for (int r = 0; r <= maxDegree; r++) {
      for (int c = 0; c <= std::min(maxOrder, r); c++) {
        //printf("(%d,%d) %.13e %.13e %.13e %.13e\n", r,c,A_real(r,c), A_imag(r,c), B_real(r,c), B_imag(r,c));
        //printf("(%d,%d) %.13e %.13e %.13e %.13e\n", r,c,A_real_ref(r,c), A_imag_ref(r,c), B_real_ref(r,c), B_imag_ref(r,c));
        assert(std::abs(A_real(r, c) - A_real_ref(r, c)) < 1e-13);
        assert(std::abs(A_imag(r, c) - A_imag_ref(r, c)) < 1e-13);
        assert(std::abs(B_real(r, c) - B_real_ref(r, c)) < 1e-13);
        assert(std::abs(B_imag(r, c) - B_imag_ref(r, c)) < 1e-13);
      }
    }

    printf("-> Test size is %dx%d : all good!\n", maxDegree, maxOrder);
  }

  return 0;
}
