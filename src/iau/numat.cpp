#include "iau.hpp"

Eigen::Matrix<double,3,3> iers2010::sofa::numat(double epsa, double dpsi,
                                  double deps) noexcept {
  //return (Eigen::AngleAxisd(epsa + deps, Eigen::Vector3d::UnitX()) *
  //        Eigen::AngleAxisd(dpsi, Eigen::Vector3d::UnitZ()) *
  //        Eigen::AngleAxisd(-epsa, Eigen::Vector3d::UnitX()))
  //    .toRotationMatrix();
  Eigen::Matrix<double,3,3> R = Eigen::Matrix<double, 3, 3>::Identity(); 
  dso::rotate<dso::RotationAxis::X>(epsa, R);
  dso::rotate<dso::RotationAxis::Z>(-dpsi, R);
  dso::rotate<dso::RotationAxis::X>(-(epsa + deps), R);
  return R;
}
