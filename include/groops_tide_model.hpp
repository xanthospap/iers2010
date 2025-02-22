#ifndef __DSO_GROOPS_TIDE_ATLAS_W_ADMITANCE_REPR_HPP__
#define __DSO_GROOPS_TIDE_ATLAS_W_ADMITANCE_REPR_HPP__

#include "eigen3/Eigen/Eigen"
#ifdef SPARSE_ADMITTANCE_MATRIX
#include "eigen3/Eigen/Sparse"
#endif

namespace dso {

class GroopsTideModel {
private:
  Eigen::Matrix<int, Eigen::Dynamic, 6> _doodson;
#ifdef SPARSE_ADMITTANCE_MATRIX
  Eigen::SparseMatrix<double> _admittance;
#else
  Eigen::MatrixXd _admittance;
#endif
public:
  GroopsTideModel() noexcept : _doodson(), _admittance() {};
#ifdef SPARSE_ADMITTANCE_MATRIX
  GroopsTideModel(const Eigen::Matrix<int, Eigen::Dynamic, 6> &d,
                  const Eigen::SparseMatrix<double> &a) noexcept
      : _doodson(d), _admittance(a) {};
  const Eigen::SparseMatrix<double> &admittance_matrix() const noexcept {
    return _admittance;
  }
#else
  GroopsTideModel(const Eigen::Matrix<int, Eigen::Dynamic, 6> &d,
                  const Eigen::MatrixXd &a) noexcept
      : _doodson(d), _admittance(a) {};
  const Eigen::MatrixXd &admittance_matrix() const noexcept {
    return _admittance;
  }
#endif
  const Eigen::Matrix<int, Eigen::Dynamic, 6> &doodson_matrix() const noexcept {
    return _doodson;
  }
  int num_waves() const noexcept { return _doodson.rows(); }
  int num_major_waves() const noexcept { return _admittance.rows(); }
}; /* class GroopsTideAtlas */

} /* namespace dso */

#endif
