#include "fcn.hpp"
#include <cstdio>
#include <vector>

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc!=2) {
    fprintf(stderr, "Usage %s [table-asc-c04.txt]\n", argv[1]);
    return 1;
  }

  /* vector of Lambert coefficients */
  std::vector<fcn::LambertCoefficients> lvec;
  
  if (parse_lambert_coefficients(argv[1], lvec)) {
    fprintf(stderr, "ERROR Failed to parse file %s\n", argv[1]);
    return 1;
  }
  printf("--> Read Lambert coefficients, size=%ld\n", lvec.size());

  for (int mjd = 51544; mjd < 53000; mjd += 10) {
    TwoPartDate t(mjd, 0e0);
    const auto res = lambert_fcn(t, lvec);
    printf("%d %.15f\n", t.imjd(), res.xcip);
  }

  return 0;
}
