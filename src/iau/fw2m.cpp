#include "iau.hpp"
#include "eigen3/Eigen/src/Geometry/RotationBase.h"

Eigen::Matrix<double,3,3> iers2010::sofa::fw2m(double gamb, double phib, double psi,
                                        double eps) noexcept {
  // r.rotz(gamb);
  // r.rotx(phib);
  // r.rotz(-psi);
  // r.rotx(-eps);
  return (Eigen::AngleAxisd(eps, Eigen::Vector3d::UnitX()) *
          Eigen::AngleAxisd(psi, Eigen::Vector3d::UnitZ()) *
          Eigen::AngleAxisd(-phib, Eigen::Vector3d::UnitX()) *
          Eigen::AngleAxisd(-gamb, Eigen::Vector3d::UnitZ()))
      .toRotationMatrix();
  /*Eigen::Matrix<double, 3, 3>::Identity();*/

  /*
   * Note that the multiplication with the identity matrix is only there to 
   * force an (implicit) cast of the Eigen::Quaternion<double> to an
   * Eigen::Matrix<double,3,3> matrix (aka, without the multiplication the 
   * expresion would return an Eigen::Quaternion<double> instance
   */
}
