#include "iau.hpp"

Eigen::Matrix<double, 3, 3> iers2010::sofa::fw2m(double gamb, double phib,
                                                 double psi,
                                                 double eps) noexcept {
  /*return (Eigen::AngleAxisd(eps, Eigen::Vector3d::UnitX()) *
          Eigen::AngleAxisd(psi, Eigen::Vector3d::UnitZ()) *
          Eigen::AngleAxisd(-phib, Eigen::Vector3d::UnitX()) *
          Eigen::AngleAxisd(-gamb, Eigen::Vector3d::UnitZ()))
      .toRotationMatrix();*/
  Eigen::Matrix<double, 3, 3> R = Eigen::Matrix<double, 3, 3>::Identity();
  dso::rotate<dso::RotationAxis::Z>(gamb, R);
  dso::rotate<dso::RotationAxis::X>(phib, R);
  dso::rotate<dso::RotationAxis::Z>(-psi, R);
  dso::rotate<dso::RotationAxis::X>(-eps, R);
  return R;
}
