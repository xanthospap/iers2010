#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
#include <utility>

using namespace dso;

/* should NOT compile! */

int main() {
  CoeffMatrix2D<MatrixStorageType::LwTriangularRowWise> mat1(4,4);
  assert(mat1.num_elements() == 10);
  assert(mat1.rows() == 4);
  assert(mat1.cols() == 4);

  [[maybe_unused]] double *col_ptr = mat1.col(0);
  return 0;
}
