#include "iau.hpp"
#include "iersc.hpp"
#include <iostream>

void apply(const iers2010::RotationMatrix3 &mat3, const double *vec, double *res) noexcept {
  res[0] = mat3.data[0][0] * vec[0] + mat3.data[0][1] * vec[1] + mat3.data[0][2]* vec[2];
  res[1] = mat3.data[1][0] * vec[0] + mat3.data[1][1] * vec[1] + mat3.data[1][2]* vec[2];
  res[2] = mat3.data[2][0] * vec[0] + mat3.data[2][1] * vec[1] + mat3.data[2][2]* vec[2];
}

int main() {
  double tt0, tt1;
  double xp, yp;
  iers2010::RotationMatrix3 rt2c;
  double xyzc[3], xyzt[3];

  while (std::cin >> tt0 >> tt1 >> xyzt[0] >> xyzt[1] >> xyzt[2]) {
    double ut1 = tt1 + 0.00789e0;
    iers2010::sofa::xy06(tt0, tt1, xp, yp);
    rt2c = iers2010::sofa::c2t06a(tt0, tt1, tt0, ut1, xp,yp);
    rt2c.transpose_inplace();
    apply(rt2c, xyzt, xyzc);
    // printf("%.9f %.3f %.3f %.3f %.3f %.3f %.3f\n", tt0+tt1, xyzt[0], xyzt[1], xyzt[2], xyzc[0], xyzc[1], xyzc[2]);
  }

  return 0;
}
