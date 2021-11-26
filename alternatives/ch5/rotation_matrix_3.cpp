#include "iau.hpp"
#include "iersc.hpp"
#include <cmath>

void iers2010::RotationMatrix3::set_identity() noexcept {
  data[0][0] = 1e0;
  data[0][1] = 0e0;
  data[0][2] = 0e0;
  data[1][0] = 0e0;
  data[1][1] = 1e0;
  data[1][2] = 0e0;
  data[2][0] = 0e0;
  data[2][1] = 0e0;
  data[2][2] = 1e0;
}

void iers2010::RotationMatrix3::rotx(double phi) noexcept {
  const double s = std::sin(phi);
  const double c = std::cos(phi);

  const double a10 = c * data[1][0] + s * data[2][0];
  const double a11 = c * data[1][1] + s * data[2][1];
  const double a12 = c * data[1][2] + s * data[2][2];
  const double a20 = -s * data[1][0] + c * data[2][0];
  const double a21 = -s * data[1][1] + c * data[2][1];
  const double a22 = -s * data[1][2] + c * data[2][2];

  data[1][0] = a10;
  data[1][1] = a11;
  data[1][2] = a12;
  data[2][0] = a20;
  data[2][1] = a21;
  data[2][2] = a22;

  return;
}

void iers2010::RotationMatrix3::roty(double phi) noexcept {
  const double s = std::sin(phi);
  const double c = std::cos(phi);

  const double a00 = c * data[0][0] - s * data[2][0];
  const double a01 = c * data[0][1] - s * data[2][1];
  const double a02 = c * data[0][2] - s * data[2][2];
  const double a20 = s * data[0][0] + c * data[2][0];
  const double a21 = s * data[0][1] + c * data[2][1];
  const double a22 = s * data[0][2] + c * data[2][2];

  data[0][0] = a00;
  data[0][1] = a01;
  data[0][2] = a02;
  data[2][0] = a20;
  data[2][1] = a21;
  data[2][2] = a22;

  return;
}

void iers2010::RotationMatrix3::rotz(double phi) noexcept {
  const double s = std::sin(phi);
  const double c = std::cos(phi);

  const double a00 = c * data[0][0] + s * data[1][0];
  const double a01 = c * data[0][1] + s * data[1][1];
  const double a02 = c * data[0][2] + s * data[1][2];
  const double a10 = -s * data[0][0] + c * data[1][0];
  const double a11 = -s * data[0][1] + c * data[1][1];
  const double a12 = -s * data[0][2] + c * data[1][2];

  data[0][0] = a00;
  data[0][1] = a01;
  data[0][2] = a02;
  data[1][0] = a10;
  data[1][1] = a11;
  data[1][2] = a12;

  return;
}

void iers2010::RotationMatrix3::mult_inplace(
    const iers2010::RotationMatrix3 &b) noexcept {
  const double c00 = data[0][0] * b.data[0][0] + data[0][1] * b.data[1][0] +
                     data[0][2] * b.data[2][0];
  const double c01 = data[0][0] * b.data[0][1] + data[0][1] * b.data[1][1] +
                     data[0][2] * b.data[2][1];
  const double c02 = data[0][0] * b.data[0][2] + data[0][1] * b.data[1][2] +
                     data[0][2] * b.data[2][2];

  const double c10 = data[1][0] * b.data[0][0] + data[1][1] * b.data[1][0] +
                     data[1][2] * b.data[2][0];
  const double c11 = data[1][0] * b.data[0][1] + data[1][1] * b.data[1][1] +
                     data[1][2] * b.data[2][1];
  const double c12 = data[1][0] * b.data[0][2] + data[1][1] * b.data[1][2] +
                     data[1][2] * b.data[2][2];

  const double c20 = data[2][0] * b.data[0][0] + data[2][1] * b.data[1][0] +
                     data[2][2] * b.data[2][0];
  const double c21 = data[2][0] * b.data[0][1] + data[2][1] * b.data[1][1] +
                     data[2][2] * b.data[2][1];
  const double c22 = data[2][0] * b.data[0][2] + data[2][1] * b.data[1][2] +
                     data[2][2] * b.data[2][2];

  data[0][0] = c00;
  data[0][1] = c01;
  data[0][2] = c02;
  data[1][0] = c10;
  data[1][1] = c11;
  data[1][2] = c12;
  data[2][0] = c20;
  data[2][1] = c21;
  data[2][2] = c22;

  return;
}

iers2010::RotationMatrix3 iers2010::RotationMatrix3::operator*(
    const iers2010::RotationMatrix3 &b) const noexcept {
  RotationMatrix3 c;

  c.data[0][0] = data[0][0] * b.data[0][0] + data[0][1] * b.data[1][0] +
                 data[0][2] * b.data[2][0];
  c.data[0][1] = data[0][0] * b.data[0][1] + data[0][1] * b.data[1][1] +
                 data[0][2] * b.data[2][1];
  c.data[0][2] = data[0][0] * b.data[0][2] + data[0][1] * b.data[1][2] +
                 data[0][2] * b.data[2][2];

  c.data[1][0] = data[1][0] * b.data[0][0] + data[1][1] * b.data[1][0] +
                 data[1][2] * b.data[2][0];
  c.data[1][1] = data[1][0] * b.data[0][1] + data[1][1] * b.data[1][1] +
                 data[1][2] * b.data[2][1];
  c.data[1][2] = data[1][0] * b.data[0][2] + data[1][1] * b.data[1][2] +
                 data[1][2] * b.data[2][2];

  c.data[2][0] = data[2][0] * b.data[0][0] + data[2][1] * b.data[1][0] +
                 data[2][2] * b.data[2][0];
  c.data[2][1] = data[2][0] * b.data[0][1] + data[2][1] * b.data[1][1] +
                 data[2][2] * b.data[2][1];
  c.data[2][2] = data[2][0] * b.data[0][2] + data[2][1] * b.data[1][2] +
                 data[2][2] * b.data[2][2];

  return c;
}

void iers2010::RotationMatrix3::transpose_inplace() noexcept {
  const double c10 = data[0][1];
  const double c20 = data[0][2];
  const double c01 = data[1][0];
  const double c21 = data[1][2];
  const double c02 = data[2][0];
  const double c12 = data[2][1];
  
  data[0][1] = c01;
  data[0][2] = c02;
  data[1][0] = c10;
  data[1][2] = c12;
  data[2][0] = c20;
  data[2][1] = c21;
}