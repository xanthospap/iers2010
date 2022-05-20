#include "matvec.hpp"
#include <cmath>

void dso::Mat3x3::set_identity() noexcept {
  data[0] = 1e0;
  data[1] = 0e0;
  data[2] = 0e0;
  data[3] = 0e0;
  data[4] = 1e0;
  data[5] = 0e0;
  data[6] = 0e0;
  data[7] = 0e0;
  data[8] = 1e0;
}

void dso::Mat3x3::rotx(double phi) noexcept {
  const double s = std::sin(phi);
  const double c = std::cos(phi);

  const double a10 = c * data[3] + s * data[6];
  const double a11 = c * data[4] + s * data[7];
  const double a12 = c * data[5] + s * data[8];
  const double a20 = -s * data[3] + c * data[6];
  const double a21 = -s * data[4] + c * data[7];
  const double a22 = -s * data[5] + c * data[8];

  data[3] = a10;
  data[4] = a11;
  data[5] = a12;
  data[6] = a20;
  data[7] = a21;
  data[8] = a22;

  return;
}

void dso::Mat3x3::roty(double phi) noexcept {
  const double s = std::sin(phi);
  const double c = std::cos(phi);

  const double a00 = c * data[0] - s * data[6];
  const double a01 = c * data[1] - s * data[7];
  const double a02 = c * data[2] - s * data[8];
  const double a20 = s * data[0] + c * data[6];
  const double a21 = s * data[1] + c * data[7];
  const double a22 = s * data[2] + c * data[8];

  data[0] = a00;
  data[1] = a01;
  data[2] = a02;
  data[6] = a20;
  data[7] = a21;
  data[8] = a22;

  return;
}

void dso::Mat3x3::rotz(double phi) noexcept {
  const double s = std::sin(phi);
  const double c = std::cos(phi);

  const double a00 = c * data[0] + s * data[3];
  const double a01 = c * data[1] + s * data[4];
  const double a02 = c * data[2] + s * data[5];
  const double a10 = -s * data[0] + c * data[3];
  const double a11 = -s * data[1] + c * data[4];
  const double a12 = -s * data[2] + c * data[5];

  data[0] = a00;
  data[1] = a01;
  data[2] = a02;
  data[3] = a10;
  data[4] = a11;
  data[5] = a12;

  return;
}

void dso::Mat3x3::mult_inplace(const dso::Mat3x3 &b) noexcept {
  const double c00 =
      data[0] * b.data[0] + data[1] * b.data[3] + data[2] * b.data[6];
  const double c01 =
      data[0] * b.data[1] + data[1] * b.data[4] + data[2] * b.data[7];
  const double c02 =
      data[0] * b.data[2] + data[1] * b.data[5] + data[2] * b.data[8];

  const double c10 =
      data[3] * b.data[0] + data[4] * b.data[3] + data[5] * b.data[6];
  const double c11 =
      data[3] * b.data[1] + data[4] * b.data[4] + data[5] * b.data[7];
  const double c12 =
      data[3] * b.data[2] + data[4] * b.data[5] + data[5] * b.data[8];

  const double c20 =
      data[6] * b.data[0] + data[7] * b.data[3] + data[8] * b.data[6];
  const double c21 =
      data[6] * b.data[1] + data[7] * b.data[4] + data[8] * b.data[7];
  const double c22 =
      data[6] * b.data[2] + data[7] * b.data[5] + data[8] * b.data[8];

  data[0] = c00;
  data[1] = c01;
  data[2] = c02;
  data[3] = c10;
  data[4] = c11;
  data[5] = c12;
  data[6] = c20;
  data[7] = c21;
  data[8] = c22;

  return;
}

dso::Mat3x3 dso::Mat3x3::operator*(const dso::Mat3x3 &b) const noexcept {
  Mat3x3 c;

  c.data[0] = data[0] * b.data[0] + data[1] * b.data[3] + data[2] * b.data[6];
  c.data[1] = data[0] * b.data[1] + data[1] * b.data[4] + data[2] * b.data[7];
  c.data[2] = data[0] * b.data[2] + data[1] * b.data[5] + data[2] * b.data[8];

  c.data[3] = data[3] * b.data[0] + data[4] * b.data[3] + data[5] * b.data[6];
  c.data[4] = data[3] * b.data[1] + data[4] * b.data[4] + data[5] * b.data[7];
  c.data[5] = data[3] * b.data[2] + data[4] * b.data[5] + data[5] * b.data[8];

  c.data[6] = data[6] * b.data[0] + data[7] * b.data[3] + data[8] * b.data[6];
  c.data[7] = data[6] * b.data[1] + data[7] * b.data[4] + data[8] * b.data[7];
  c.data[8] = data[6] * b.data[2] + data[7] * b.data[5] + data[8] * b.data[8];

  return c;
}

dso::Mat3x3& dso::Mat3x3::transpose_inplace() noexcept {
  const double c10 = data[1];
  const double c20 = data[2];
  const double c01 = data[3];
  const double c21 = data[5];
  const double c02 = data[6];
  const double c12 = data[7];

  data[1] = c01;
  data[2] = c02;
  data[3] = c10;
  data[5] = c12;
  data[6] = c20;
  data[7] = c21;

  return *this;
}

dso::Mat3x3 dso::Mat3x3::transpose() noexcept {
  const double c10 = data[1];
  const double c20 = data[2];
  const double c01 = data[3];
  const double c21 = data[5];
  const double c02 = data[6];
  const double c12 = data[7];

  return dso::Mat3x3{{data[0], c01, c02, c10, data[4], c12, c20, c21, data[8]}};
}
