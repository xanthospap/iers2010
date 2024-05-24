#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
using namespace dso;

int main() {
  [[maybe_unused]]CoeffMatrix2D<MatrixStorageType::ColumnWise> mat1(500,400);
  [[maybe_unused]]CoeffMatrix2D<MatrixStorageType::ColumnWise> mat2(500,400);
  [[maybe_unused]]CoeffMatrix2D<MatrixStorageType::ColumnWise> mat3(5,4);

  /* not allowed! */
  mat3 += 1e0 * mat1.reduced_view(321, 299) + 2e0 * mat2.reduced_view(321, 299);

  return 0;
};
