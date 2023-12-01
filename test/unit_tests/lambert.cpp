#include "fcn.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage %s [table-asc-c04.txt]\n", argv[1]);
    return 1;
  }

  /* vector of Lambert coefficients */
  std::vector<fcn::LambertCoefficients> lvec;

  /* parse all data from input file */
  if (parse_lambert_coefficients(argv[1], lvec)) {
    fprintf(stderr, "ERROR Failed to parse file %s\n", argv[1]);
    return 1;
  }

  TwoPartDate t;
  for (double mjd = 51544e0; mjd < 53000e0; mjd += 10.03) {
    /* make the date (split parts) */
    {
      double imjd, secday;
      secday = std::modf(mjd, &imjd) * dso::SEC_PER_DAY;
      t = TwoPartDate((int)imjd, secday);
    }
    const auto res = lambert_fcn(t, lvec);
    printf("%20.5f %18.12f %18.12f %18.12f %18.12f\n", t.as_mjd(), res.xcip,
           res.ycip, res.sxcip, res.sycip);
  }

  return 0;
}
