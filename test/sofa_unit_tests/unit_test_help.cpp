#include "unit_test_help.hpp"
#include <limits>

// Seed with a real random value, if available
std::random_device r;

// random engine
std::default_random_engine e(r());

bool approx_equal(double a, double b, double epsilon) {
  return (a == b) || (std::abs(a - b) < epsilon);
}

bool approx_equal(double a, double b) {
  return a == b || std::abs(a - b) < std::abs(std::min(a, b)) *
                                         std::numeric_limits<double>::epsilon();
}

bool approx_equal(double a, double b, const char *error_msg) {
  bool equal = (a == b) ||
               (std::abs(a - b) < std::abs(std::min(a, b)) *
                                      std::numeric_limits<double>::epsilon());
  if (!equal) {
    fprintf(stderr, "ERROR. values: %.12e != %.12e, diff=%.6e\n", a, b, a - b);
    fprintf(stderr, "%s\n", error_msg);
  }

  return equal;
}

bool approx_equal(double a, double b, double epsilon, const char *error_msg) {
  bool equal = (a == b) || (std::abs(a - b) < epsilon);

  if (!equal) {
    fprintf(stderr, "ERROR. values: %.12e != %.12e, diff=%.6e, Epsilon=%.12e\n",
            a, b, a - b, epsilon);
    fprintf(stderr, "%s\n", error_msg);
  }

  return equal;
}

bool approx_equal(const dso::TwoPartDate &d1, const dso::TwoPartDate &d2) {
  const double a = d1.big() + d1.small();
  const double b = d2.big() + d2.small();
  return a == b || std::abs(a - b) < std::abs(std::min(a, b)) *
                                         std::numeric_limits<double>::epsilon();
}

bool approx_equal(const Eigen::Matrix<double, 3, 3> &m1, const double m2[3][3],
                  bool print_on_fail) noexcept {
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
  if (equal != 9 && print_on_fail) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++)
        printf("%+20.12f ", m1(i, j));
      printf("\t");
      for (int j = 0; j < 3; j++)
        printf("%+20.12f ", m2[i][j]);
      printf("\t");
      for (int j = 0; j < 3; j++)
        printf("%+5e ", m1(i, j) - m2[i][j]);
      printf("\n");
    }
  }
  return equal == 9;
}

bool approx_equal(const Eigen::Matrix<double, 3, 3> &m1, const double m2[3][3],
                  double epsilon, bool print_on_fail) noexcept {
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
  if (equal != 9 && print_on_fail) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++)
        printf("%+20.12f ", m1(i, j));
      printf("\t");
      for (int j = 0; j < 3; j++)
        printf("%+20.12f ", m2[i][j]);
      printf("\t");
      for (int j = 0; j < 3; j++)
        printf("%+5e ", m1(i, j) - m2[i][j]);
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
  return dso::TwoPartDate(d.big(), d.small() + uni_sec(e) / 86400e0)
      .normalized();
}

double random_angle(double min, double max) noexcept {
  std::uniform_real_distribution<double> uni_rad(min, max);
  return uni_rad(e);
}

/* returns angle in [rad] */
double rotation_matrix_diff(const Eigen::Matrix<double, 3, 3> &m1,
                            const double m2[3][3]) {
  Eigen::Matrix<double, 3, 3> m22;
  m22 << m2[0][0], m2[0][1], m2[0][2], m2[1][0], m2[1][1], m2[1][2], m2[2][0],
      m2[2][1], m2[2][2];

  /* see
   * https://math.stackexchange.com/questions/2113634/comparing-two-rotation-matrices
   * and https://www.cs.cmu.edu/~cga/dynopt/readings/Rmetric.pdf
   * Ra <- m1
   * Rb <- m2
   */
  const auto Rab = m1.transpose() * m22;

  /* Tr(Rab) = 1 + 2 cos(θ)
   * θ = arccos( (Tr(Rab) - 1) / 2 )
   */
  if (std::abs((Rab.trace() - 1e0) / 2.) > 1e0) {
    if (std::abs(Rab.trace()-3e0) > 1e-16) {
    fprintf(stderr,
            "#ERROR@%s Invalid arccos argument: %.6e outside [-1,1] range "
            "(trace=%.3e)\n",
            __func__, (Rab.trace() - 1e0) / 2., Rab.trace());
    return std::numeric_limits<double>::max();
    } else {
      ;
    }
  }
  const double theta = std::acos((Rab.trace() - 1e0) / 2.);
  return theta;
}
