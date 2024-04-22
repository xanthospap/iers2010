#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
#include <utility>

using namespace dso;

void dprint(const CoeffMatrix2D<MatrixStorageType::ColumnWise> &mat) {
  for (int r=0; r<mat.rows(); r++) {
    for (int c=0; c<mat.cols(); c++) {
      printf("%+6.2f ", mat(r,c));
    }
    printf("\n");
  }
  return;
}

int main() {
  CoeffMatrix2D<MatrixStorageType::ColumnWise> mat1(5,4);
  assert(mat1.num_elements() == 20);
  assert(mat1.rows() == 5);
  assert(mat1.cols() == 4);

  /* fill in */
  int k=0;
  for (int i=0; i<mat1.rows(); i++) {
    for (int j=0; j<mat1.cols(); j++) {
      mat1(i,j) = (double)(++k);
    }
  }
  dprint(mat1);

  /* resize copying data that are there */
  mat1.cresize(10,5);
  dprint(mat1);
  {
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

  /* get columns */
  double *cptr = mat1.column(0);
  assert(*cptr == mat1(0,0));
  assert(*(cptr+1) == mat1(1,0));
  assert(*(cptr+2) == mat1(2,0));
  assert(*(cptr+3) == mat1(3,0));
  assert(*(cptr+4) == mat1(4,0));
  cptr = mat1.column(1);
  assert(*cptr == mat1(0,1));
  assert(*(cptr+1) == mat1(1,1));
  assert(*(cptr+2) == mat1(2,1));
  assert(*(cptr+3) == mat1(3,1));
  assert(*(cptr+4) == mat1(4,1));
  cptr = mat1.column(2);
  assert(*cptr == mat1(0,2));
  assert(*(cptr+1) == mat1(1,2));
  assert(*(cptr+2) == mat1(2,2));
  assert(*(cptr+3) == mat1(3,2));
  assert(*(cptr+4) == mat1(4,2));
  cptr = mat1.column(3);
  assert(*cptr == mat1(0,3));
  assert(*(cptr+1) == mat1(1,3));
  assert(*(cptr+2) == mat1(2,3));
  assert(*(cptr+3) == mat1(3,3));
  assert(*(cptr+4) == mat1(4,3));
  }

  /* resize copying data that are there */
  mat1.cresize(3,4);
  dprint(mat1);
  {
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

  /* get columns */
  double *cptr = mat1.column(0);
  assert(*cptr == mat1(0,0));
  assert(*(cptr+1) == mat1(1,0));
  assert(*(cptr+2) == mat1(2,0));
  cptr = mat1.column(1);
  assert(*cptr == mat1(0,1));
  assert(*(cptr+1) == mat1(1,1));
  assert(*(cptr+2) == mat1(2,1));
  cptr = mat1.column(2);
  assert(*cptr == mat1(0,2));
  assert(*(cptr+1) == mat1(1,2));
  assert(*(cptr+2) == mat1(2,2));
  cptr = mat1.column(3);
  assert(*cptr == mat1(0,3));
  assert(*(cptr+1) == mat1(1,3));
  assert(*(cptr+2) == mat1(2,3));
  }
  
  return 0;
}
