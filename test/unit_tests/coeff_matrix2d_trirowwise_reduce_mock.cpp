#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
#include <utility>

using namespace dso;

int main() {
  CoeffMatrix2D<MatrixStorageType::LwTriangularRowWise> mat1(5,5);
  CoeffMatrix2D<MatrixStorageType::LwTriangularRowWise> mat2(2,2);
  
  /* fill in */
  //int k=0;
  //for (int i=0; i<mat1.rows(); i++) {
  //  for (int j=0; j<=i; j++) {
  //    mat1(i,j) = (double)(++k);
  //    mat2(i,j) = (double)(++k) * .1;
  //  }
  //}

  CoeffMatrix2D<MatrixStorageType::LwTriangularRowWise> mat3 =
      1e0 * mat1.reduced_view(3, 3) + 2e0 * mat2.reduced_view(3, 3);
  assert(mat3.rows() == 3);
  assert(mat3.cols() == 3);
  //for (int i=0; i<3; i++) {
  //  for (int j=0; j<=i; j++) {
  //    assert(mat3(i, j) == 1e0 * mat1(i, j) + 2e0 * mat2(i, j));
  //  }
  //}

  return 0;
}
