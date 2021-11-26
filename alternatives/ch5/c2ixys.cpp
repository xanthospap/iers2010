#include "iau.hpp"
#include "iersc.hpp"
#include <cmath>

void iers2010::iau::c2ixys(double x, double y, double s, iers2010::RotationMatrix3 &r3mat) noexcept {
  
  // Obtain the spherical angles E and d
  const double r2 = x*x +y*y;
  const double e = (r2>0e0) ? std::atan2(y,x) : 0e0;
  const double d = std::atan(std::sqrt(r2 / (1e0-r2)));

    // initialize rotation matrix to identity
    r3mat.set_identity();

    // R3 (−E) x R2 (−d) x R3 (E) x R3 (s),
    // X = sin d cos E,
    // Y = sin d sin E, 
    // Z = cos d, 
    r3mat.rotz(e);
    r3mat.roty(d);
    r3mat.rotz(-(e+s));

    return;
}