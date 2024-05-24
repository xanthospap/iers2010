#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
#include <utility>

using namespace dso;

int main() {
  CoeffMatrix2D<MatrixStorageType::LwTriangularRowWise> mat1(5,5);
  CoeffMatrix2D<MatrixStorageType::LwTriangularRowWise> mat2(2,2);
  
  /* this line should fail! */
  CoeffMatrix2D<MatrixStorageType::LwTriangularRowWise> mat3 =
      1e0 * mat1.reduced_view(3, 3) + 2e0 * mat2.reduced_view(3, 3);
  
  assert(mat3.rows() == 3);
  assert(mat3.cols() == 3);

  return 0;
}
