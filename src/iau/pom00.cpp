#include "iau.hpp"

Eigen::Matrix<double, 3, 3> iers2010::sofa::pom00(double xp, double yp,
                                                  double sp) noexcept {
  
  //return
  //    (Eigen::AngleAxisd(-yp, -Eigen::Vector3d::UnitX()) *
  //     Eigen::AngleAxisd(-xp, -Eigen::Vector3d::UnitY()) *
  //     Eigen::AngleAxisd(sp, -Eigen::Vector3d::UnitZ()))
  //        .toRotationMatrix();
  
  /* apply three rotations ... RPOM(t) = R3(−sp) x R2(xp) x R1(yp) */
  Eigen::Matrix<double, 3, 3> R2 = Eigen::Matrix<double, 3, 3>::Identity();
  dso::rotate<dso::RotationAxis::Z>(sp,R2);
  dso::rotate<dso::RotationAxis::Y>(-xp,R2);
  dso::rotate<dso::RotationAxis::X>(-yp,R2);
  return R2;
}
