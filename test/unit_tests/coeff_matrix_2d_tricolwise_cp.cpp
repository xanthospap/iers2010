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

void dprint(const CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &mat) {
  for (int r = 0; r < mat.rows(); r++) {
    for (int c = 0; c <= r; c++) {
      printf("%+6.2f ", mat(r, c));
    }
    printf("\n");
  }
  return;
}

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
    mat2.cresize(2, 2);
    assert(mat2(0, 0) == 1);
    assert(mat2(1, 0) == 5);
    assert(mat2(1, 1) == 6);

    /* get columns */
    double *cptr = mat2.column(0);
    assert(*cptr == mat2(0, 0));
    assert(*(cptr + 1) == mat2(1, 0));
    cptr = mat2.column(1);
    assert(*cptr == mat2(1, 1));
  }

  {
    auto mat2(mat1);
    mat2.cresize(8, 8);
    assert(mat2(0, 0) == 1);
    assert(mat2(1, 0) == 5);
    assert(mat2(2, 0) == 9);
    assert(mat2(3, 0) == 13);
    assert(mat2(1, 1) == 6);
    assert(mat2(2, 1) == 10);
    assert(mat2(3, 1) == 14);
    assert(mat2(2, 2) == 11);
    assert(mat2(3, 2) == 15);
    assert(mat2(3, 3) == 16);

    /* get columns */
    double *cptr = mat2.column(0);
    assert(*cptr == mat2(0, 0));
    assert(*(cptr + 1) == mat2(1, 0));
    assert(*(cptr + 2) == mat2(2, 0));
    assert(*(cptr + 3) == mat2(3, 0));
    cptr = mat2.column(1);
    assert(*cptr == mat2(1, 1));
    assert(*(cptr + 1) == mat2(2, 1));
    assert(*(cptr + 2) == mat2(3, 1));
    cptr = mat2.column(2);
    assert(*cptr == mat2(2, 2));
    assert(*(cptr + 1) == mat2(3, 2));
    cptr = mat2.column(3);
    assert(*cptr == mat2(3, 3));
  }
  
  {
    auto mat2(mat1);
    mat2.cresize(1, 1);
    assert(mat2(0, 0) == 1);

    /* get columns */
    double *cptr = mat2.column(0);
    assert(*cptr == mat2(0, 0));
  }

  /* all done */
  return 0;
}
