#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 10000;

const char *funcs[] = {"pn00"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);
const char *args[] = {"rb", "rp", "rbp", "rn", "rbpn", "epsa"};
const int num_args = sizeof(args) / sizeof(args[0]);

int main() {
  int fails[6] = {0, 0, 0, 0, 0, 0};
  int error = 0;
  double max_error[6] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
  double rb[3][3], rp[3][3], rbp[3][3], rn[3][3], rbpn[3][3];
  Eigen::Matrix<double, 3, 3> rbm, rpm, rbpm, rnm, rbpnm;
  double epsas, epsam;

  printf("Function         #Tests #Fails #Maxerror[masec]  Status\n");
  printf("---------------------------------------------------------------\n");

  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const double dpsi = random_angle(-iers2010::DPI / 5, iers2010::DPI / 5);
    const double deps = random_angle(-iers2010::DPI / 5, iers2010::DPI / 5);
    pn00(tt, dpsi, deps, epsam, rbm, rpm, rbpm, rnm, rbpnm);
    iauPn00(tt.big() + dso::mjd0_jd, tt.small(), dpsi, deps, &epsas, rb, rp,
            rbp, rn, rbpn);
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
        max_error[1] = theta;
      }
    }
    /* rbp */
    if (!approx_equal(rbpm, rbp)) {
      const double theta = rotation_matrix_diff(rbpm, rbp);
      ++fails[2];
      if (std::abs(theta) > std::abs(max_error[2])) {
        max_error[2] = theta;
      }
    }
    /* rn */
    if (!approx_equal(rnm, rn)) {
      const double theta = rotation_matrix_diff(rnm, rn);
      ++fails[3];
      if (std::abs(theta) > std::abs(max_error[3])) {
        max_error[3] = theta;
      }
    }
    /* rbpn */
    if (!approx_equal(rbpnm, rbpn)) {
      const double theta = rotation_matrix_diff(rbpnm, rbpn);
      ++fails[4];
      if (std::abs(theta) > std::abs(max_error[4])) {
        max_error[4] = theta;
      }
    }
    /* epsa */
    if (!approx_equal(epsam, epsas)) {
      ++fails[5];
      if (std::abs(epsam - epsas) > std::abs(max_error[5])) {
        max_error[5] = epsas - epsam;
      }
    }
  }

  error = 0;
  for (int j = 0; j < 6; j++) {
    printf("%8s %7s %6d %6d %+.1e %.7s %s\n", funcs[0], args[j], NUM_TESTS,
           fails[j], dso::rad2sec(max_error[j])*1e3,
           (fails[j] == 0) ? "OK" : "FAILED",
           (j == 5) ? "Angle" : "RotMatrix");
    error += fails[j];
  }

  return error;
}
