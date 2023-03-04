#ifndef __IERS2010_UTEST_IAU_HPP__
#define __IERS2010_UTEST_IAU_HPP__

#include "eigen3/Eigen/Eigen"
#include "iau.hpp"
#include "sofa.h"
#include <limits>
#include <random>

inline bool approx_equal(double a, double b) {
  return a == b || std::abs(a - b) < std::abs(std::min(a, b)) *
                                         std::numeric_limits<double>::epsilon();
}

bool approx_equal(double a, double b, const char *errormsg);

inline bool approx_equal(double a, double b, double epsilon) {
  return (a == b) || (std::abs(a - b) < epsilon);
}

bool approx_equal(double a, double b, double epsilon, const char *errormsg);

bool approx_equal(const dso::TwoPartDate &d1, const dso::TwoPartDate &d2);

bool approx_equal(const Eigen::Matrix<double, 3, 3> &m1,
                  const double m2[3][3], bool print_on_fail=false) noexcept;

bool approx_equal(const Eigen::Matrix<double, 3, 3> &m1, const double m2[3][3],
                  double epsilon, bool print_on_fail = false) noexcept;

dso::TwoPartDate random_mjd();

dso::TwoPartDate add_random_seconds(const dso::TwoPartDate &d,
                                    double min_sec = -60e0 * 10,
                                    double max_sec = 60e0 * 10);

double random_angle(double min = -dso::D2PI, double max = dso::D2PI) noexcept;

#endif
