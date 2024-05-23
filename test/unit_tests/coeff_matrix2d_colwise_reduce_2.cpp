#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
#include <utility>

using namespace dso;

void printmat(const CoeffMatrix2D<MatrixStorageType::ColumnWise> &mat) {
  for (int i=0; i<mat.rows(); i++) {
    for (int j=0; j<mat.cols(); j++) {
      printf("%6.2f ", mat(i,j));
    }
    printf("\n");
  }
  return;
}

int main() {
  CoeffMatrix2D<MatrixStorageType::ColumnWise> mat1(500,400);
  CoeffMatrix2D<MatrixStorageType::ColumnWise> mat2(500,400);
  
  /* fill in */
  int k=0;
  for (int i=0; i<mat1.rows(); i++) {
    for (int j=0; j<mat1.cols(); j++) {
      mat1(i,j) = (double)(++k);
      mat2(i,j) = (double)(++k) * .1;
    }
  }

  CoeffMatrix2D<MatrixStorageType::ColumnWise> mat3 =
      1e0 * mat1.reduced_view(321, 299) + 2e0 * mat2.reduced_view(321, 299);
  
  /* validation */
  assert(mat3.rows() == 321);
  assert(mat3.cols() == 299);
  for (int i=0; i<321; i++) {
    for (int j=0; j<299; j++) {
      assert(mat3(i, j) == 1e0 * mat1(i, j) + 2e0 * mat2(i, j));
    }
  }

  mat3 += 1e0 * mat1.reduced_view(321, 299) + 2e0 * mat2.reduced_view(321, 299);

  /* validation */
  assert(mat3.rows() == 321);
  assert(mat3.cols() == 299);
  for (int i=0; i<321; i++) {
    for (int j=0; j<299; j++) {
      assert(mat3(i, j) == 2e0*(1e0 * mat1(i, j) + 2e0 * mat2(i, j)));
    }
  }

  return 0;
}
