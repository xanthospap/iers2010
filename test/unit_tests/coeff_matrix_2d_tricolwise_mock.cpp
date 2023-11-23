#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
#include <utility>

using namespace dso;

/* This should not compile !
 */

int main() {
  CoeffMatrix2D<MatrixStorageType::ColumnWise> mat1(5,4);
  assert(mat1.num_elements() == 20);
  
  [[maybe_unused]] double *row_ptr = mat1.row(0);
  
  return 0;
}