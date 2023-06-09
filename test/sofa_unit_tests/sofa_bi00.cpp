#include "geodesy/units.hpp"
#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010::sofa;

constexpr const int NUM_TESTS = 10000;

const char *funcs[] = {"bi00"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);
const char *args[] = {"dpsibi", "depsbi", "dra"};
const int num_args = sizeof(args) / sizeof(args[0]);

int main(int argc, [[maybe_unused]] char *argv[]) {
  if (argc > 1) {
    fprintf(stderr, "Ignoring command line arguments!\n");
  }
  int func_it = 0;
  int fails[3] = {0,0,0};
  int error = 0;
  double max_error[3] = {0e0,0e0,0e0};
  double am[3], as[3];

  printf("Function         #Tests #Fails #Maxerror[sec]    Status\n");
  printf("---------------------------------------------------------------\n");

  for (int i = 0; i < NUM_TESTS; i++) {
    bi00(am[0], am[1], am[2]);
    iauBi00(&as[0], &as[1], &as[2]);
    for (int j = 0; j < 3; j++) {
      if (!approx_equal(am[j], as[j])) {
        ++fails[j];
        if (std::abs(am[j] - as[j]) > max_error[j]) {
          max_error[j] = std::abs(am[j] - as[j]);
        }
      }
    }
  }

  for (int j = 0; j < 3; j++) {
    printf("%8s %7s %6d %6d %+.9e %s\n", funcs[func_it], args[j], NUM_TESTS,
           fails[j], dso::rad2sec(max_error[j]), (fails[j] == 0) ? "OK" : "FAILED");
  }

  error = 0;
  for (int j = 0; j < 3; j++) error += fails[j];

  return error;
}
