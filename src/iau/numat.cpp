#include "iau.hpp"
#include "eigen3/Eigen/src/Geometry/RotationBase.h"

Eigen::Matrix<double,3,3> iers2010::sofa::numat(double epsa, double dpsi,
                                  double deps) noexcept {
  //Eigen::Matrix<double,3,3> rmatn;
  //rmatn.rotx(epsa);
  //rmatn.rotz(-dpsi);
  //rmatn.rotx(-(epsa + deps));
  //return rmatn;
  return (Eigen::AngleAxisd(epsa + deps, Eigen::Vector3d::UnitX()) *
          Eigen::AngleAxisd(dpsi, Eigen::Vector3d::UnitZ()) *
          Eigen::AngleAxisd(-epsa, Eigen::Vector3d::UnitX()))
      .toRotationMatrix();
}
