#ifndef __DSO_COSTG_BENCHMARK_UTILS_HPP__
#define __DSO_COSTG_BENCHMARK_UTILS_HPP__

#include <vector>
#include "eigen3/Eigen/Eigen"
#include "datetime/calendar.hpp"

namespace costg {

struct BmOrbit {
  dso::MjdEpoch epoch;
  Eigen::Matrix<double, 3, 1> xyz;
  Eigen::Matrix<double, 3, 1> vxyz;
  Eigen::Matrix<double, 3, 1> axyz;
  BmOrbit(const dso::MjdEpoch &t, const double *rva) : epoch(t) {
    xyz << rva[0], rva[1], rva[2];
    vxyz << rva[3], rva[4], rva[5];
    axyz << rva[6], rva[7], rva[8];
  }
};

struct BmAcceleration {
  dso::MjdEpoch epoch;
  Eigen::Matrix<double, 3, 1> axyz;
  BmAcceleration(const dso::MjdEpoch &t, double ax, double ay, double az)
      : epoch(t) {
    axyz << ax, ay, az;
  }
};

struct BmRotaryMatrix {
  dso::MjdEpoch epoch;
  Eigen::Matrix<double, 3, 3> R;
  BmRotaryMatrix(const dso::MjdEpoch &t, const double *data)
      : epoch(t) {
    R << data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8];
  }
};

struct BmEops {
  dso::MjdEpoch epoch;
  double xp, yp, sp, dUT1, LOD, X, Y, s;
  BmEops(const dso::MjdEpoch &t, const double *data)
      : epoch(t), xp(data[0]), yp(data[1]), sp(data[2]), dUT1(data[3]),
        LOD(data[4]), X(data[5]), Y(data[6]), s(data[7]){};
};

  const char *skipws(const char *line) noexcept;

  std::vector<BmAcceleration> parse_acceleration(const char *fn);

  std::vector<BmOrbit> parse_orbit(const char *fn);

  std::vector<BmRotaryMatrix> parse_rotary(const char *fn);

  std::vector<BmEops> parse_eops(const char *fn);

  const char *basename(const char *fn);

} /* namespace costg */

#endif
