#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
#include <utility>

using namespace dso;

/*
 * Matrix1 (rows=5, cols=4)
 *  1  2  3  4
 *  5  6  7  8
 *  9 10 11 12
 * 13 14 15 16
 * 17 18 19 20
 */

int main() {
  CoeffMatrix2D<MatrixStorageType::RowWise> mat1(5,4);
  assert(mat1.num_elements() == 20);
  
  /* fill in */
  int k=0;
  for (int i=0; i<mat1.rows(); i++) {
    for (int j=0; j<mat1.cols(); j++) {
      mat1(i,j) = (double)(++k);
    }
  }

  assert(mat1(0,0) == 1);
  assert(mat1(0,1) == 2);
  assert(mat1(0,2) == 3);
  assert(mat1(0,3) == 4);
  assert(mat1(1,0) == 5);
  assert(mat1(1,1) == 6);
  assert(mat1(1,2) == 7);
  assert(mat1(1,3) == 8);
  assert(mat1(2,0) == 9);
  assert(mat1(2,1) == 10);
  assert(mat1(2,2) == 11);
  assert(mat1(2,3) == 12);
  assert(mat1(3,0) == 13);
  assert(mat1(3,1) == 14);
  assert(mat1(3,2) == 15);
  assert(mat1(3,3) == 16);
  assert(mat1(4,0) == 17);
  assert(mat1(4,1) == 18);
  assert(mat1(4,2) == 19);
  assert(mat1(4,3) == 20);

  /* get rows */
  double *rptr = mat1.row(0);
  assert(*rptr == mat1(0,0));
  assert(*(rptr+1) == mat1(0,1));
  assert(*(rptr+2) == mat1(0,2));
  assert(*(rptr+3) == mat1(0,3));
  rptr = mat1.row(1);
  assert(*rptr == mat1(1,0));
  assert(*(rptr+1) == mat1(1,1));
  assert(*(rptr+2) == mat1(1,2));
  assert(*(rptr+3) == mat1(1,3));
  rptr = mat1.row(2);
  assert(*rptr == mat1(2,0));
  assert(*(rptr+1) == mat1(2,1));
  assert(*(rptr+2) == mat1(2,2));
  assert(*(rptr+3) == mat1(2,3));
  rptr = mat1.row(3);
  assert(*rptr == mat1(3,0));
  assert(*(rptr+1) == mat1(3,1));
  assert(*(rptr+2) == mat1(3,2));
  assert(*(rptr+3) == mat1(3,3));
  rptr = mat1.row(4);
  assert(*rptr == mat1(4,0));
  assert(*(rptr+1) == mat1(4,1));
  assert(*(rptr+2) == mat1(4,2));
  assert(*(rptr+3) == mat1(4,3));

  /* scale */
  mat1.multiply(2e0);
  assert(mat1(0,0) == 2e0*1);
  assert(mat1(0,1) == 2e0*2);
  assert(mat1(0,2) == 2e0*3);
  assert(mat1(0,3) == 2e0*4);
  assert(mat1(1,0) == 2e0*5);
  assert(mat1(1,1) == 2e0*6);
  assert(mat1(1,2) == 2e0*7);
  assert(mat1(1,3) == 2e0*8);
  assert(mat1(2,0) == 2e0*9);
  assert(mat1(2,1) == 2e0*10);
  assert(mat1(2,2) == 2e0*11);
  assert(mat1(2,3) == 2e0*12);
  assert(mat1(3,0) == 2e0*13);
  assert(mat1(3,1) == 2e0*14);
  assert(mat1(3,2) == 2e0*15);
  assert(mat1(3,3) == 2e0*16);
  assert(mat1(4,0) == 2e0*17);
  assert(mat1(4,1) == 2e0*18);
  assert(mat1(4,2) == 2e0*19);
  assert(mat1(4,3) == 2e0*20);

  /* copy */
  auto mat2(mat1);
  assert(mat2(0,0) == 2e0*1);
  assert(mat2(0,1) == 2e0*2);
  assert(mat2(0,2) == 2e0*3);
  assert(mat2(0,3) == 2e0*4);
  assert(mat2(1,0) == 2e0*5);
  assert(mat2(1,1) == 2e0*6);
  assert(mat2(1,2) == 2e0*7);
  assert(mat2(1,3) == 2e0*8);
  assert(mat2(2,0) == 2e0*9);
  assert(mat2(2,1) == 2e0*10);
  assert(mat2(2,2) == 2e0*11);
  assert(mat2(2,3) == 2e0*12);
  assert(mat2(3,0) == 2e0*13);
  assert(mat2(3,1) == 2e0*14);
  assert(mat2(3,2) == 2e0*15);
  assert(mat2(3,3) == 2e0*16);
  assert(mat2(4,0) == 2e0*17);
  assert(mat2(4,1) == 2e0*18);
  assert(mat2(4,2) == 2e0*19);
  assert(mat2(4,3) == 2e0*20);
  
  auto mat3=mat2;
  assert(mat3(0,0) == 2e0*1);
  assert(mat3(0,1) == 2e0*2);
  assert(mat3(0,2) == 2e0*3);
  assert(mat3(0,3) == 2e0*4);
  assert(mat3(1,0) == 2e0*5);
  assert(mat3(1,1) == 2e0*6);
  assert(mat3(1,2) == 2e0*7);
  assert(mat3(1,3) == 2e0*8);
  assert(mat3(2,0) == 2e0*9);
  assert(mat3(2,1) == 2e0*10);
  assert(mat3(2,2) == 2e0*11);
  assert(mat3(2,3) == 2e0*12);
  assert(mat3(3,0) == 2e0*13);
  assert(mat3(3,1) == 2e0*14);
  assert(mat3(3,2) == 2e0*15);
  assert(mat3(3,3) == 2e0*16);
  assert(mat3(4,0) == 2e0*17);
  assert(mat3(4,1) == 2e0*18);
  assert(mat3(4,2) == 2e0*19);
  assert(mat3(4,3) == 2e0*20);
  
  auto mat4(std::move(mat3));
  assert(mat4(0,0) == 2e0*1);
  assert(mat4(0,1) == 2e0*2);
  assert(mat4(0,2) == 2e0*3);
  assert(mat4(0,3) == 2e0*4);
  assert(mat4(1,0) == 2e0*5);
  assert(mat4(1,1) == 2e0*6);
  assert(mat4(1,2) == 2e0*7);
  assert(mat4(1,3) == 2e0*8);
  assert(mat4(2,0) == 2e0*9);
  assert(mat4(2,1) == 2e0*10);
  assert(mat4(2,2) == 2e0*11);
  assert(mat4(2,3) == 2e0*12);
  assert(mat4(3,0) == 2e0*13);
  assert(mat4(3,1) == 2e0*14);
  assert(mat4(3,2) == 2e0*15);
  assert(mat4(3,3) == 2e0*16);
  assert(mat4(4,0) == 2e0*17);
  assert(mat4(4,1) == 2e0*18);
  assert(mat4(4,2) == 2e0*19);
  assert(mat4(4,3) == 2e0*20);
  
  assert(mat3.num_elements() == 0);

  /* set */
  mat4.fill_with(-5e0);
  for (int i=0; i<mat1.rows(); i++) {
    for (int j=0; j<mat1.cols(); j++) {
      assert(mat4(i,j) == -5e0);
    }
  }

  /* all done */
  return 0;
}
