#include "iau.hpp"
#include "eigen3/Eigen/src/Geometry/RotationBase.h"

Eigen::Matrix<double,3,3> iers2010::sofa::pom00(double xp, double yp, double sp) noexcept {
  // apply three rotations ... W(t) = R3(−sp) x R2(xp) x R1(yp),
  //rpom.rotz(sp);
  //rpom.roty(-xp);
  //rpom.rotx(-yp);
  //return rpom;
  return (Eigen::AngleAxisd(yp, Eigen::Vector3d::UnitX()) *
          Eigen::AngleAxisd(xp, Eigen::Vector3d::UnitY()) *
          Eigen::AngleAxisd(-sp, Eigen::Vector3d::UnitZ()))
      .toRotationMatrix();
}
