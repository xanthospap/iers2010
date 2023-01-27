#include "unit_test_help.hpp"

// Seed with a real random value, if available
std::random_device r;

// random engine
std::default_random_engine e(r());

bool approx_equal(const dso::TwoPartDate &d1, const dso::TwoPartDate &d2) {
  const double a = d1._big + d1._small;
  const double b = d2._big + d2._small;
  return a == b || std::abs(a - b) < std::abs(std::min(a, b)) *
                                         std::numeric_limits<double>::epsilon();
}

bool approx_equal(const Eigen::Matrix<double, 3, 3> &m1,
                  const double m2[3][3], bool print_on_fail) noexcept {
  int equal = 0;
  equal += approx_equal(m1(0, 0), m2[0][0]);
  equal += approx_equal(m1(1, 0), m2[1][0]);
  equal += approx_equal(m1(2, 0), m2[2][0]);
  equal += approx_equal(m1(0, 1), m2[0][1]);
  equal += approx_equal(m1(1, 1), m2[1][1]);
  equal += approx_equal(m1(2, 1), m2[2][1]);
  equal += approx_equal(m1(0, 2), m2[0][2]);
  equal += approx_equal(m1(1, 2), m2[1][2]);
  equal += approx_equal(m1(2, 2), m2[2][2]);
  if (equal!=9 && print_on_fail) {
    for (int i=0; i<3; i++) {
      for (int j=0;j<3;j++)
        printf("%+20.12f ", m1(i,j));
      printf("\t");
      for (int j=0;j<3;j++)
        printf("%+20.12f ", m2[i][j]);
      printf("\t");
      for (int j=0;j<3;j++)
        printf("%+5e ", m1(i,j)-m2[i][j]);
      printf("\n");
  }
  }
  return equal == 9;
}

bool approx_equal(const Eigen::Matrix<double, 3, 3> &m1,
                  const double m2[3][3], double epsilon, bool print_on_fail) noexcept {
  int equal = 0;
  equal += approx_equal(m1(0, 0), m2[0][0], epsilon);
  equal += approx_equal(m1(1, 0), m2[1][0], epsilon);
  equal += approx_equal(m1(2, 0), m2[2][0], epsilon);
  equal += approx_equal(m1(0, 1), m2[0][1], epsilon);
  equal += approx_equal(m1(1, 1), m2[1][1], epsilon);
  equal += approx_equal(m1(2, 1), m2[2][1], epsilon);
  equal += approx_equal(m1(0, 2), m2[0][2], epsilon);
  equal += approx_equal(m1(1, 2), m2[1][2], epsilon);
  equal += approx_equal(m1(2, 2), m2[2][2], epsilon);
  if (equal!=9 && print_on_fail) {
    for (int i=0; i<3; i++) {
      for (int j=0;j<3;j++)
        printf("%+20.12f ", m1(i,j));
      printf("\t");
      for (int j=0;j<3;j++)
        printf("%+20.12f ", m2[i][j]);
      printf("\t");
      for (int j=0;j<3;j++)
        printf("%+5e ", m1(i,j)-m2[i][j]);
      printf("\n");
  }
  }
  return equal == 9;
}

dso::TwoPartDate random_mjd() {
  constexpr const double min_mjd = 29629e0; // 01/01/1940
  constexpr const double max_mjd = 69807e0; // 01/01/2050
  std::uniform_real_distribution<double> uni_mjd(min_mjd, max_mjd);
  std::uniform_real_distribution<double> uni_fmjd(0e0, 1e0);
  return dso::TwoPartDate(uni_mjd(e), uni_fmjd(e)).normalized();
}

dso::TwoPartDate add_random_seconds(const dso::TwoPartDate &d, double min_sec,
                                    double max_sec) {
  std::uniform_real_distribution<double> uni_sec(min_sec, max_sec);
  return dso::TwoPartDate(d._big, d._small + uni_sec(e) / 86400e0).normalized();
}

double random_angle(double min, double max) noexcept {
  std::uniform_real_distribution<double> uni_rad(min, max);
  return uni_rad(e);
}
