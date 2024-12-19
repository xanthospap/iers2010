#include "gravity.hpp"

Eigen::Matrix<double, 3, 1>
dso::point_mass_acceleration(const Eigen::Matrix<double, 3, 1> &r,
                             const Eigen::Matrix<double, 3, 1> &rcb,
                             double GMcb) noexcept {

  const double Rcb3 = rcb.squaredNorm() * rcb.norm();
  const auto s = r - rcb;
  const double S3 = s.squaredNorm() * s.norm();

  return (s / S3 + rcb / Rcb3) * (-GMcb);
}

Eigen::Matrix<double, 3, 1>
dso::point_mass_acceleration(const Eigen::Matrix<double, 3, 1> &r,
                             const Eigen::Matrix<double, 3, 1> &rcb,
                             double GMcb,
                             Eigen::Matrix<double, 3, 3> &J) noexcept {

  const double Rcb3 = rcb.squaredNorm() * rcb.norm();
  const auto s = r - rcb;
  const double S3 = s.squaredNorm() * s.norm();
  const double S5 = S3 * s.squaredNorm();

  /* jacobian */
  J = (-GMcb) * (Eigen::Matrix<double, 3, 3>::Identity() / S3 -
                 3e0 * s * s.transpose() / S5);

  /* acceleration vector */
  return (s / S3 + rcb / Rcb3) * (-GMcb);
}
