#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
#include <utility>

using namespace dso;

/* Should NOT compile !
 */

int main() {
  CoeffMatrix2D<MatrixStorageType::RowWise> mat1(5,4);
  assert(mat1.num_elements() == 20);

  [[maybe_unused]] double *col_ptr = mat1.col(0);

  return 0;
}
