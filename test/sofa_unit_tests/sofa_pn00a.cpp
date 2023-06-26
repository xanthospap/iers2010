#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"pn00a"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);
const char *args[] = {"dpsi", "deps", "epsa", "rb", "rp", "rbp", "rn", "rbpn"};
const char *type[] = {"Angle",       "Angle",       "Angle",
                      "RotMatrix", "RotMatrix", "RotMatrix",
                      "RotMatrix", "RotMatrix"};
const int num_args = sizeof(args) / sizeof(args[0]);

int main() {
  int fails[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int error = 0;
  double max_error[8] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
  double rb[3][3], rp[3][3], rbp[3][3], rn[3][3], rbpn[3][3];
  Eigen::Matrix<double, 3, 3> rbm, rpm, rbpm, rnm, rbpnm;
  double dpsis, dpsim, depsm, depss, epsas, epsam;

  printf("Function         #Tests #Fails #Maxerror[masec]  Status  Type\n");
  printf("---------------------------------------------------------------\n");

  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    pn00a(tt, dpsim, depsm, epsam, rbm, rpm, rbpm, rnm, rbpnm);
    iauPn00a(tt.big() + dso::mjd0_jd, tt.small(), &dpsis, &depss, &epsas, rb,
             rp, rbp, rn, rbpn);
    /* dpsi */
    if (!approx_equal(dpsim, dpsis)) {
      ++fails[0];
      if (std::abs(dpsim - dpsis) > std::abs(max_error[0])) {
        max_error[0] = dpsis - dpsim;
      }
    }
    /* deps */
    if (!approx_equal(depsm, depss)) {
      ++fails[1];
      if (std::abs(depsm - depss) > std::abs(max_error[1])) {
        max_error[1] = depss - depsm;
      }
    }
    /* epsa */
    if (!approx_equal(epsam, epsas)) {
      ++fails[2];
      if (std::abs(epsam - epsas) > std::abs(max_error[2])) {
        max_error[2] = epsas - epsam;
      }
    }
    /* rb */
    if (!approx_equal(rbm, rb)) {
      ++fails[3];
      const double theta = rotation_matrix_diff(rbm, rb);
      if (std::abs(theta) > std::abs(max_error[3])) {
        max_error[3] = theta;
      }
    }
    /* rp */
    if (!approx_equal(rpm, rp)) {
      ++fails[4];
      const double theta = rotation_matrix_diff(rpm, rp);
      if (std::abs(theta) > std::abs(max_error[4])) {
        max_error[4] = theta;
      }
    }
    /* rbp */
    if (!approx_equal(rbpm, rbp)) {
      ++fails[5];
      const double theta = rotation_matrix_diff(rbpm, rbp);
      if (std::abs(theta) > std::abs(max_error[5])) {
        max_error[5] = theta;
      }
    }
    /* rn */
    if (!approx_equal(rnm, rn)) {
      ++fails[6];
      const double theta = rotation_matrix_diff(rnm, rn);
      if (std::abs(theta) > std::abs(max_error[6])) {
        max_error[6] = theta;
      }
    }
    /* rbpn */
    if (!approx_equal(rbpnm, rbpn)) {
      ++fails[7];
      const double theta = rotation_matrix_diff(rbpnm, rbpn);
      if (std::abs(theta) > std::abs(max_error[7])) {
        max_error[7] = theta;
      }
    }
  }

  error = 0;
  for (int j = 0; j < 8; j++) {
    printf("%8s %7s %6d %6d %+.1e %.7s %s\n", funcs[0], args[j], NUM_TESTS,
           fails[j], dso::rad2sec(max_error[j])*1e3,
           (fails[j] == 0) ? "OK" : "FAILED", type[j]);
    error += fails[j];
  }

  return error;
}
