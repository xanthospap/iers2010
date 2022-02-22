#ifndef __MATVEC_3X3_SIMPLE_HPP__
#define __MATVEC_3X3_SIMPLE_HPP__

#include <cmath>
#include <cstring>

namespace dso {

struct Vector3 {
  double data[3] = {0e0, 0e0, 0e0};
  constexpr const double &x() const noexcept {return data[0];}
  constexpr const double &y() const noexcept {return data[1];}
  constexpr const double &z() const noexcept {return data[2];}
  constexpr double &x() noexcept {return data[0];}
  constexpr double &y() noexcept {return data[1];}
  constexpr double &z() noexcept {return data[2];}

  static constexpr Vector3 to_vec3(const double *vec) noexcept {
    return Vector3{{vec[0], vec[1], vec[2]}};
  }

  constexpr double norm_squared() const noexcept {
    return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
  }
  constexpr double norm() const noexcept {
    return std::sqrt(this->norm_squared());
  }
  constexpr double dot_product(const Vector3 &v) const noexcept {
    return data[0] * v.data[0] + data[1] * v.data[1] + data[2] * v.data[2];
  }
  constexpr Vector3 cross_product(const Vector3 &v) const noexcept {
    const double s1 = data[1] * v.data[2] - data[2] * v.data[1];
    const double s2 = data[2] * v.data[0] - data[0] * v.data[2];
    const double s3 = data[0] * v.data[1] - data[1] * v.data[0];
    return Vector3{{s1, s2, s3}};
  }
  constexpr void cross_product(const Vector3 &v, Vector3 &out) const noexcept {
    out.data[0] = data[1] * v.data[2] - data[2] * v.data[1];
    out.data[1] = data[2] * v.data[0] - data[0] * v.data[2];
    out.data[2] = data[0] * v.data[1] - data[1] * v.data[0];
  }
  constexpr Vector3 normalize() const noexcept {
    const double norm = this->norm();
    return Vector3{{data[0] / norm, data[1] / norm, data[2] / norm}};
  }
  constexpr void normalize() noexcept {
    const double norm = this->norm();
    data[0] /= norm;
    data[1] /= norm;
    data[2] /= norm;
  }
  constexpr Vector3 operator-(const Vector3 &v) const noexcept {
    return {{data[0]-v.data[0], data[1]-v.data[1], data[2]-v.data[2]}};
  }
  constexpr Vector3 operator+(const Vector3 &v) const noexcept {
    return {{data[0]+v.data[0], data[1]+v.data[1], data[2]+v.data[2]}};
  }
  constexpr Vector3 operator/(double scalar) const noexcept {
    return {{data[0]/scalar, data[1]/scalar, data[2]/scalar}};
  }
  constexpr Vector3 operator*(double scalar) const noexcept {
    return {{data[0]*scalar, data[1]*scalar, data[2]*scalar}};
  }

};

struct Mat3x3 {
  double data[9] = {1e0, 0e0, 0e0, 0e0, 1e0, 0e0, 0e0, 0e0, 1e0};

  constexpr double &operator()(int i, int j) noexcept {
    return data[i * 3 + j];
  }
  constexpr const double &operator()(int i, int j) const noexcept {
    return data[i * 3 + j];
  }

  /// @brief Multiply two 3x3 matrices (aka this * b)
  Mat3x3 operator*(const Mat3x3 &b) const noexcept;

  /// @brief Multiply a vector by a matrix.
  Vector3 operator*(const Vector3 &vec) const noexcept {
    double v0 =
        data[0] * vec.data[0] + data[1] * vec.data[1] + data[2] * vec.data[2];
    double v1 =
        data[3] * vec.data[0] + data[4] * vec.data[1] + data[5] * vec.data[2];
    double v2 =
        data[6] * vec.data[0] + data[7] * vec.data[1] + data[8] * vec.data[2];
    return Vector3{{v0, v1, v2}};
  }

  /// @brief Multiply two 3x3 matrices (aka this * b) and store result in
  ///        this instance
  void mult_inplace(const Mat3x3 &b) noexcept;

  /// @brief Transpose a 3x3 matric (in place)
  void transpose_inplace() noexcept;

  /// @brief Transpose a 3x3 matric
  Mat3x3 transpose() noexcept;

  /// @brief Set to identity matrix
  void set_identity() noexcept;

  /// @brief Rotate an r-matrix about the x-axis.
  /// This will actually perform the operation R = Rx * R, with Rx =
  ///  1        0            0
  ///  0   + cos(phi)   + sin(phi)
  ///  0   - sin(phi)   + cos(phi)
  /// @param[in] angle (radians)
  void rotx(double angle) noexcept;

  /// @brief Rotate an r-matrix about the y-axis.
  /// This will actually perform the operation R = Ry * R, with Rx =
  ///  + cos(phi)     0      - sin(phi)
  ///       0           1           0
  ///  + sin(phi)     0      + cos(phi)
  /// @param[in] angle (radians)
  void roty(double) noexcept;

  /// @brief Rotate an r-matrix about the z-axis.
  /// This will actually perform the operation R = Rz * R, with Rx =
  ///  + cos(psi)   + sin(psi)     0
  ///  - sin(psi)   + cos(psi)     0
  ///       0            0         1
  /// @param[in] angle (radians)
  void rotz(double) noexcept;
}; // Mat3x3
} // dso

#endif
