#include "iau.hpp"
#include "iersc.hpp"
#include <iostream>

int main() {
  double tt0, tt1, x, y, z;
  double xp, yp;
  iers2010::RotationMatrix3 rt2c;

  while (std::cin >> tt0 >> tt1 >> x >>y >> z) {
    double ut1 = tt1 + 0.00789e0;
    iers2010::sofa::xy06(tt0, tt1, xp, yp);
    rt2c = iers2010::sofa::c2t06a(tt0, tt1, tt0, ut1, xp,yp);
  }

  return 0;
}
