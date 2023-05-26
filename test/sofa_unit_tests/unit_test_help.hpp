#ifndef __IERS2010_UTEST_IAU_HPP__
#define __IERS2010_UTEST_IAU_HPP__

#include "eigen3/Eigen/Eigen"
#include "iau.hpp"
#include "sofa.h"
#include <limits>
#include <random>

bool approx_equal(double a, double b);

bool approx_equal(double a, double b, const char *errormsg);

bool approx_equal(double a, double b, double epsilon);

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

double rotation_matrix_diff(const Eigen::Matrix<double, 3, 3> &m1, const double m2[3][3]);

inline Eigen::Matrix<double,3,3>
rxr2mat(const double a[3][3]) {
  Eigen::Matrix<double,3,3> m;
  m(0,0) = a[0][0];
  m(0,1) = a[0][1];
  m(0,2) = a[0][2];
  m(1,0) = a[1][0];
  m(1,1) = a[1][1];
  m(1,2) = a[1][2];
  m(2,0) = a[2][0];
  m(2,1) = a[2][1];
  m(2,2) = a[2][2];
  return m;
}

#endif
