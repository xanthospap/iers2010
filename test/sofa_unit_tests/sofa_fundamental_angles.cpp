#include "iau.hpp"
#include "sofa.h"
#include "unit_test_help.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

using namespace iers2010;

constexpr const int NUM_TESTS = 100000;

const char *funcs[] = {"fapa03", "fane03", "faur03", "fasa03", "faju03",
                       "fama03", "fae03",  "fave03", "fame03", "faom03",
                       "fad03",  "faf03",  "falp03", "fal03"};
const int num_funs = sizeof(funcs) / sizeof(funcs[0]);

int main() {
  int fails[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double max_error[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int error = 0;

  printf("Function         #Tests #Fails #Maxerror[sec]    Status  Type\n");
  printf("---------------------------------------------------------------\n");

  double am, as;
  int j = 0;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fapa03(tt);
    as = iauFapa03(jc);
    if (!approx_equal(am, as)) {
      ++fails[j];
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fane03(tt);
    as = iauFane03(jc);
    if (!approx_equal(am, as)) {
      ++fails[j];
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = faur03(tt);
    as = iauFaur03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fasa03(tt);
    as = iauFasa03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = faju03(tt);
    as = iauFaju03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fama03(tt);
    as = iauFama03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fae03(tt);
    as = iauFae03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fave03(tt);
    as = iauFave03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fame03(tt);
    as = iauFame03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = faom03(tt);
    as = iauFaom03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fad03(tt);
    as = iauFad03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = faf03(tt);
    as = iauFaf03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = falp03(tt);
    as = iauFalp03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  ++j;
  for (int i = 0; i < NUM_TESTS; i++) {
    const auto tt = random_mjd();
    const auto jc = tt.jcenturies_sinceJ2000();
    am = fal03(tt);
    as = iauFal03(jc);
    if (!approx_equal(am, as)) {
      if (std::abs(am - as) > max_error[j])
        max_error[j] = std::abs(am - as);
      ++fails[j];
    }
  }
  printf("%8s %7s %6d %6d %+.9e %.7s %s\n", funcs[j], funcs[j], NUM_TESTS,
         fails[j], dso::rad2sec(max_error[j]),
         (fails[j] == 0) ? "OK" : "FAILED", "Angle");
  if (fails[j])
    ++error;

  return error;
}
