#ifndef __DSO_ROTATION_MATRICES_HPP__
#define __DSO_ROTATION_MATRICES_HPP__

#include "eigen3/Eigen/Eigen"

namespace dso {

enum class RotationAxis : char {X,Y,Z};

namespace {
void
rotateX(double angle, Eigen::Matrix<double, 3, 3> &R) noexcept {
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  const double a10 =   c*R(1,0) + s*R(2,0);
  const double a11 =   c*R(1,1) + s*R(2,1);
  const double a12 =   c*R(1,2) + s*R(2,2);
  const double a20 = - s*R(1,0) + c*R(2,0);
  const double a21 = - s*R(1,1) + c*R(2,1);
  const double a22 = - s*R(1,2) + c*R(2,2);

   R(1,0) = a10;
   R(1,1) = a11;
   R(1,2) = a12;
   R(2,0) = a20;
   R(2,1) = a21;
   R(2,2) = a22;
};

void
rotateY(double angle, Eigen::Matrix<double, 3, 3> &R) noexcept {
   const double s = std::sin(angle);
   const double c = std::cos(angle);

   const double a00 = c*R(0,0) - s*R(2,0);
   const double a01 = c*R(0,1) - s*R(2,1);
   const double a02 = c*R(0,2) - s*R(2,2);
   const double a20 = s*R(0,0) + c*R(2,0);
   const double a21 = s*R(0,1) + c*R(2,1);
   const double a22 = s*R(0,2) + c*R(2,2);

   R(0,0) = a00;
   R(0,1) = a01;
   R(0,2) = a02;
   R(2,0) = a20;
   R(2,1) = a21;
   R(2,2) = a22;
};

void
rotateZ(double angle, Eigen::Matrix<double, 3, 3> &R) noexcept {
   const double s = std::sin(angle);
   const double c = std::cos(angle);

   const double a00 =   c*R(0,0) + s*R(1,0);
   const double a01 =   c*R(0,1) + s*R(1,1);
   const double a02 =   c*R(0,2) + s*R(1,2);
   const double a10 = - s*R(0,0) + c*R(1,0);
   const double a11 = - s*R(0,1) + c*R(1,1);
   const double a12 = - s*R(0,2) + c*R(1,2);

   R(0,0) = a00;
   R(0,1) = a01;
   R(0,2) = a02;
   R(1,0) = a10;
   R(1,1) = a11;
   R(1,2) = a12;
};
} /* unnamed namespace */

template <RotationAxis T>
void rotate(double angle, Eigen::Matrix<double, 3, 3> &R) noexcept {
   if constexpr (T == RotationAxis::X)
     rotateX(angle, R);
   else if (T == RotationAxis::Y)
     rotateY(angle, R);
   else
     rotateZ(angle, R);
};

template <RotationAxis T>
[[nodiscard]] Eigen::Matrix<double, 3, 3>
rotate(double angle, const Eigen::Matrix<double, 3, 3> &Rin) noexcept {
   Eigen::Matrix<double, 3, 3> Rout(Rin);
   if constexpr (T == RotationAxis::X)
     rotateX(angle, Rout);
   else if (T == RotationAxis::Y)
     rotateY(angle, Rout);
   else
     rotateZ(angle, Rout);

   return Rout;
};

template <RotationAxis T>
[[nodiscard]] Eigen::Matrix<double, 3, 1>
rotate(double angle, const Eigen::Matrix<double, 3, 1> &r) noexcept {
   Eigen::Matrix<double, 3, 3> R = Eigen::Matrix<double, 3, 3>::Identity();
   if constexpr (T == RotationAxis::X)
     rotateX(angle, R);
   else if (T == RotationAxis::Y)
     rotateY(angle, R);
   else
     rotateZ(angle, R);

   return R*r;
};

}/* namespace dso */

#endif
