#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
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
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> mat1(4,4);
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

  assert(mat1(0,0) == 1);
  assert(mat1(1,0) == 5);
  assert(mat1(2,0) == 9);
  assert(mat1(3,0) == 13);
  assert(mat1(1,1) == 6);
  assert(mat1(2,1) == 10);
  assert(mat1(3,1) == 14);
  assert(mat1(2,2) == 11);
  assert(mat1(3,2) == 15);
  assert(mat1(3,3) == 16);

  /* get columns */
  double *cptr = mat1.column(0);
  assert(*cptr == mat1(0,0));
  assert(*(cptr+1) == mat1(1,0));
  assert(*(cptr+2) == mat1(2,0));
  assert(*(cptr+3) == mat1(3,0));
  cptr = mat1.column(1);
  assert(*cptr == mat1(1,1));
  assert(*(cptr+1) == mat1(2,1));
  assert(*(cptr+2) == mat1(3,1));
  cptr = mat1.column(2);
  assert(*cptr == mat1(2,2));
  assert(*(cptr+1) == mat1(3,2));
  cptr = mat1.column(3);
  assert(*cptr == mat1(3,3));

  /* scale */
  mat1.multiply(2e0);
  assert(mat1(0,0) == 2e0*1);
  assert(mat1(1,0) == 2e0*5);
  assert(mat1(1,1) == 2e0*6);
  assert(mat1(2,0) == 2e0*9);
  assert(mat1(2,1) == 2e0*10);
  assert(mat1(2,2) == 2e0*11);
  assert(mat1(3,0) == 2e0*13);
  assert(mat1(3,1) == 2e0*14);
  assert(mat1(3,2) == 2e0*15);
  assert(mat1(3,3) == 2e0*16);

  /* copy */
  auto mat2(mat1);
  assert((mat1.rows() == mat2.rows()) && (mat1.cols() == mat2.cols()));
  for (int i=0; i<mat2.rows(); i++) {
    for (int j=0; j<mat2.cols(); j++) {
      if (j<=i) {
        assert(mat2(i,j)==mat1(i,j));
      }
    }
  }
  
  auto mat3=mat2;
  assert((mat3.rows() == mat2.rows()) && (mat3.cols() == mat2.cols()));
  for (int i=0; i<mat3.rows(); i++) {
    for (int j=0; j<mat3.cols(); j++) {
      if (j<=i) {
        assert(mat3(i,j)==mat2(i,j));
      }
    }
  }
  
  auto mat4(std::move(mat3));
  assert((mat4.rows() == mat2.rows()) && (mat4.cols() == mat2.cols()));
  for (int i=0; i<mat4.rows(); i++) {
    for (int j=0; j<mat4.cols(); j++) {
      if (j<=i) {
        assert(mat4(i,j)==mat2(i,j));
      }
    }
  }
  assert(mat3.num_elements() == 0);

  /* set */
  mat4.fill_with(-5e0);
  for (int i=0; i<mat1.rows(); i++) {
    for (int j=0; j<mat1.cols(); j++) {
      if (j<=i) assert(mat4(i,j) == -5e0);
    }
  }

  /* all done */
  return 0;
}
