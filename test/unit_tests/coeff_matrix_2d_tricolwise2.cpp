#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
#include <utility>

using namespace dso;

/*
 * Matrix1 (rows=4, cols=4)
 *  1                        3
 *  5  6                     4  5
 *  9 10 11                  3  4  5
 * 13 14 15 16               4  6  8  10
 */

int main() {
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> mat1(4,4);
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> mat2(4,4);
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> mat3(4,4);
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> mat4(4,4);

  /* fill in */
  int k=0;
  for (int i=0; i<mat1.rows(); i++) {
    for (int j=0; j<mat1.cols(); j++) {
      ++k;
      if (j<=i) {
        mat1(i,j) = (double)(k);
      }
    }
  }

  mat2 = mat1;
  mat3(0,0) = 3;
  mat3(1,0) = 4;
  mat3(1,1) = 5;
  mat3(2,0) = 3;
  mat3(2,1) = 4;
  mat3(2,2) = 5;
  mat3(3,0) = 4;
  mat3(3,1) = 6;
  mat3(3,2) = 8;
  mat3(3,3) = 10;

  mat4 = mat1 + .1 * mat2 + .2 * mat3;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (j <= i) {
        assert(mat4(i,j) == mat1(i,j) + .1 * mat2(i,j) + .2 * mat3(i,j));
      }
    }
  }
  
  mat4 += mat1 + .1 * mat2 + .2 * mat3;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (j <= i) {
        assert(mat4(i,j) == 2e0*(mat1(i,j) + .1 * mat2(i,j) + .2 * mat3(i,j)));
      }
    }
  }

  /* all done */
  return 0;
}
