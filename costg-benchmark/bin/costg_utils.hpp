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
  BmOrbit(const dso::MjdEpoch &t, double *rva) : epoch(t) {
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

const char *skipws(const char *line) noexcept;

std::vector<BmAcceleration> parse_acceleration(const char *fn);

std::vector<BmOrbit> parse_orbit(const char *fn);

} /* namespace costg */

#endif
