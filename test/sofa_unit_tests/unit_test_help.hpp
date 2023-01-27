#include "iau.hpp"
#include "sofa.h"
#include <limits>
#include <random>
#include "eigen3/Eigen/Eigen"

// Seed with a real random value, if available
std::random_device r;

// random engine
std::default_random_engine e(r());

inline bool approx_equal(double a, double b) {
  return a == b || std::abs(a - b) < std::abs(std::min(a, b)) *
                                         std::numeric_limits<double>::epsilon();
}

inline bool approx_equal(const dso::TwoPartDate &d1,
                         const dso::TwoPartDate &d2) {
  const double a = d1._big + d1._small;
  const double b = d2._big + d2._small;
  return a == b || std::abs(a - b) < std::abs(std::min(a, b)) *
                                         std::numeric_limits<double>::epsilon();
}

inline bool approx_equal(const Eigen::Matrix<double,3,3> &m1, const double m2[3][3]) noexcept {
    int equal = 0;
    equal += approx_equal(m1(0,0), m2[0][0]);
    equal += approx_equal(m1(1,0), m2[1][0]);
    equal += approx_equal(m1(2,0), m2[2][0]);
    equal += approx_equal(m1(0,1), m2[0][1]);
    equal += approx_equal(m1(1,1), m2[1][1]);
    equal += approx_equal(m1(2,1), m2[2][1]);
    equal += approx_equal(m1(0,2), m2[0][2]);
    equal += approx_equal(m1(1,2), m2[1][2]);
    equal += approx_equal(m1(2,2), m2[2][2]);
    return equal == 9;
}

inline dso::TwoPartDate random_mjd(const std::default_random_engine &re=e) {
  constexpr const double min_mjd = 29629e0; // 01/01/1940
  constexpr const double max_mjd = 69807e0; // 01/01/2050
  std::uniform_real_distribution<double> uni_mjd(min_mjd, max_mjd);
  std::uniform_real_distribution<double> uni_fmjd(0e0, 1e0);
  return dso::TwoPartDate(uni_mjd(re), uni_fmjd(re)).normalized();
}

inline dso::TwoPartDate add_random_seconds(const dso::TwoPartDate &d,
                          double min_sec = -60e0 * 10,
                          double max_sec = 60e0 * 10) {
  std::uniform_real_distribution<double> uni_sec(min_sec, max_sec);
  return dso::TwoPartDate(d._big, d._small + uni_sec(e) / 86400e0).normalized();
}

inline double random_angle(double min=-dso::D2PI, double max=dso::D2PI) noexcept {
    std::uniform_real_distribution<double> uni_rad(min,max);
    return uni_rad(e);
}