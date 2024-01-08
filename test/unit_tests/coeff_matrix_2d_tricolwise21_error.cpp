#include "coeff_matrix_2d.hpp"
#include <cstdio>
#include <cassert>
#include <utility>

using namespace dso;

int main() {
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> mat1(4,4);
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> mat2(4,4);
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> mat3(4,4);
  CoeffMatrix2D<MatrixStorageType::LwTriangularRowWise> mat4(4,4);

  /* this should NOT be allowed to compile cause mat4 is of different type */
  mat4 = mat1 + .1 * mat2 + .2 * mat3;
  
  /* all done */
  return 0;
}
