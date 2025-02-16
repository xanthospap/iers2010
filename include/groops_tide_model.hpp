#ifndef __DSO_GROOPS_TIDE_ATLAS_W_ADMITANCE_REPR_HPP__
#define __DSO_GROOPS_TIDE_ATLAS_W_ADMITANCE_REPR_HPP__

#include "eigen3/Eigen/Eigen"

namespace dso {

class GroopsTideModel {
private:
  Eigen::Matrix<int, Eigen::Dynamic, 6> _doodson;
  Eigen::MatrixXd _admittance;
public:
  GroopsTideModel() noexcept : _doodson(), _admittance() {};
  GroopsTideModel(const Eigen::Matrix<int, Eigen::Dynamic, 6> &d, const Eigen::MatrixXd &a) noexcept : _doodson(d), _admittance(a) {};
  const Eigen::Matrix<int, Eigen::Dynamic, 6> &doodson_matrix() const noexcept { return _doodson; }
  const Eigen::MatrixXd &admittance_matrix() const noexcept { return _admittance;}
  int num_waves() const noexcept { return _doodson.rows(); }
  int num_major_waves() const noexcept { return _admittance.rows(); }
}; /* class GroopsTideAtlas */

} /* namespace dso */

#endif
