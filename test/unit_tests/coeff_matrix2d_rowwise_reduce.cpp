#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
#include <utility>

using namespace dso;

void printmat(const CoeffMatrix2D<MatrixStorageType::RowWise> &mat) {
  for (int i=0; i<mat.rows(); i++) {
    for (int j=0; j<mat.cols(); j++) {
      printf("%6.2f ", mat(i,j));
    }
    printf("\n");
  }
  return;
}

int main() {
  CoeffMatrix2D<MatrixStorageType::RowWise> mat1(5,4);
  CoeffMatrix2D<MatrixStorageType::RowWise> mat2(5,4);
  
  /* fill in */
  int k=0;
  for (int i=0; i<mat1.rows(); i++) {
    for (int j=0; j<mat1.cols(); j++) {
      mat1(i,j) = (double)(++k);
      mat2(i,j) = (double)(++k) * .1;
    }
  }

  CoeffMatrix2D<MatrixStorageType::RowWise> mat3 =
      1e0 * mat1.reduced_view(3, 2) + 2e0 * mat2.reduced_view(3, 2);
  
  assert(mat3.rows() == 3);
  assert(mat3.cols() == 2);
  for (int i=0; i<3; i++) {
    for (int j=0; j<2; j++) {
      assert(mat3(i, j) == 1e0 * mat1(i, j) + 2e0 * mat2(i, j));
    }
  }

  mat3 += 1e0 * mat1.reduced_view(3, 2) + 2e0 * mat2.reduced_view(3, 2);
  
  assert(mat3.rows() == 3);
  assert(mat3.cols() == 2);
  for (int i=0; i<3; i++) {
    for (int j=0; j<2; j++) {
      assert(mat3(i, j) == 2e0*(1e0 * mat1(i, j) + 2e0 * mat2(i, j)));
    }
  }

  return 0;
}
