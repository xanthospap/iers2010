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
  CoeffMatrix2D<MatrixStorageType::LwTriangularRowWise> mat1(4,4);
  assert(mat1.num_elements() == 10);
  assert(mat1.rows() == 4);
  assert(mat1.cols() == 4);
  
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

  {
    auto mat2(mat1);
    mat2.cresize(5,5);
  assert(mat1(0,0) == 1);
  assert(mat1(1,0) == 5);
  assert(mat1(1,1) == 6);
  assert(mat1(2,0) == 9);
  assert(mat1(2,1) == 10);
  assert(mat1(2,2) == 11);
  assert(mat1(3,0) == 13);
  assert(mat1(3,1) == 14);
  assert(mat1(3,2) == 15);
  assert(mat1(3,3) == 16);

  /* get columns */
  double *cptr = mat1.row(0);
  assert(*cptr == mat1(0,0));
  cptr = mat1.row(1);
  assert(*cptr == mat1(1,0));
  assert(cptr[1] == mat1(1,1));
  cptr = mat1.row(2);
  assert(*cptr == mat1(2,0));
  assert(cptr[1] == mat1(2,1));
  assert(cptr[2] == mat1(2,2));
  cptr = mat1.row(3);
  assert(*cptr == mat1(3,0));
  assert(cptr[1] == mat1(3,1));
  assert(cptr[2] == mat1(3,2));
  assert(cptr[3] == mat1(3,3));
  }
  
  {
    auto mat2(mat1);
    mat2.cresize(50,50);
  assert(mat1(0,0) == 1);
  assert(mat1(1,0) == 5);
  assert(mat1(1,1) == 6);
  assert(mat1(2,0) == 9);
  assert(mat1(2,1) == 10);
  assert(mat1(2,2) == 11);
  assert(mat1(3,0) == 13);
  assert(mat1(3,1) == 14);
  assert(mat1(3,2) == 15);
  assert(mat1(3,3) == 16);

  /* get columns */
  double *cptr = mat1.row(0);
  assert(*cptr == mat1(0,0));
  cptr = mat1.row(1);
  assert(*cptr == mat1(1,0));
  assert(cptr[1] == mat1(1,1));
  cptr = mat1.row(2);
  assert(*cptr == mat1(2,0));
  assert(cptr[1] == mat1(2,1));
  assert(cptr[2] == mat1(2,2));
  cptr = mat1.row(3);
  assert(*cptr == mat1(3,0));
  assert(cptr[1] == mat1(3,1));
  assert(cptr[2] == mat1(3,2));
  assert(cptr[3] == mat1(3,3));
  }
  
  {
    auto mat2(mat1);
    mat2.cresize(3,3);
  assert(mat1(0,0) == 1);
  assert(mat1(1,0) == 5);
  assert(mat1(1,1) == 6);
  assert(mat1(2,0) == 9);
  assert(mat1(2,1) == 10);
  assert(mat1(2,2) == 11);

  /* get columns */
  double *cptr = mat1.row(0);
  assert(*cptr == mat1(0,0));
  cptr = mat1.row(1);
  assert(*cptr == mat1(1,0));
  assert(cptr[1] == mat1(1,1));
  cptr = mat1.row(2);
  assert(*cptr == mat1(2,0));
  assert(cptr[1] == mat1(2,1));
  assert(cptr[2] == mat1(2,2));
  }

  return 0;
}
