#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"bp00"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);
const char *args[] = {"rb", "rp", "rbp"};
const char *type[] = {"RotMatrix", "RotMatrix", "RotMatrix"};
const int num_args = sizeof(args) / sizeof(args[0]);

int main(int argc, [[maybe_unused]] char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  int fails[] = {0, 0, 0};
  int error = 0;
  double max_error[] = {0e0, 0e0, 0e0};
  double rb[3][3], rp[3][3], rbp[3][3];
  Eigen::Matrix<double, 3, 3> rbm, rpm, rbpm;

  printf("Function         #Tests #Fails #Maxerror[masec]  Status  Type\n");
  printf("---------------------------------------------------------------\n");

  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    bp00(tt, rbm, rpm, rbpm);
    iauBp00(tt.big() + dso::mjd0_jd, tt.small(), rb, rp, rbp);
    /* rb */
    if (!approx_equal(rbm, rb)) {
      const double theta = rotation_matrix_diff(rbm, rb);
      ++fails[0];
      if (std::abs(theta) > std::abs(max_error[0])) {
        max_error[0] = theta;
      }
    }
    /* rp */
    if (!approx_equal(rpm, rp)) {
      const double theta = rotation_matrix_diff(rpm, rp);
      ++fails[1];
      if (std::abs(theta) > std::abs(max_error[1])) {
        // printf("-->New theta Angle=%.15e, t=%.20e / %.3f\n",
        // dso::rad2sec(theta), tt.small(), tt.big());
        max_error[1] = theta;
        // for (int k=0; k<3; k++) {
        //   for (int j=0;j<3;j++) {
        //     printf("%+.20e ", rpm(k,j));
        //   }
        //   printf("\n");
        // }
        // for (int k=0; k<3; k++) {
        //   for (int j=0;j<3;j++) {
        //     printf("%+.20e ", rp[k][j]);
        //   }
        //   printf("\n");
        // }
      }
    }
    /* rbp */
    if (!approx_equal(rbpm, rbp)) {
      const double theta = rotation_matrix_diff(rbpm, rbp);
      ++fails[2];
      if (std::abs(theta) > max_error[2]) {
        max_error[2] = std::abs(theta);
      }
    }
  }

  for (int j = 0; j < 3; j++) {
    printf("%8s %7s %6d %6d %+.1e %.7s %s\n", funcs[0], args[j], NUM_TESTS,
           fails[j], dso::rad2sec(max_error[j])*1e3,
           (fails[j] == 0) ? "OK" : "FAILED", type[j]);
  }

  error = 0;
  for (int j = 0; j < 3; j++)
    error += fails[j];

  return error;
}
