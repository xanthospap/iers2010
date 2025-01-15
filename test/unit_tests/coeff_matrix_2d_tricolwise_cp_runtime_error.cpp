#include "coeff_matrix_2d.hpp"
#include <cassert>
#include <cstdio>
#include <utility>

using namespace dso;

/*
 * Matrix1 (rows=4, cols=4)
 *  1
 *  5  6
 *  9 10 11
 * 13 14 15 16
 */

int main() {
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> mat1(4, 4);
  assert(mat1.num_elements() == 10);
  assert(mat1.rows() == 4);
  assert(mat1.cols() == 4);

  /* fill in */
  int k = 0;
  for (int i = 0; i < mat1.rows(); i++) {
    for (int j = 0; j < mat1.cols(); j++) {
      ++k;
      if (j <= i) {
        mat1(i, j) = (double)(k);
      }
    }
  }

  {
    auto mat2(mat1);
    /* cannot resize if rows !=  cols 
     * should call abort at runtime
     */
    mat2.cresize(4, 3);
    assert(mat2(0, 0) == 1);

    /* get columns */
    double *cptr = mat2.column(0);
    assert(*cptr == mat2(0, 0));
  }

  /* all done */
  return 0;
}
