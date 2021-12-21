#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

int main() {

  // testing orthoeop
  std::cout << "----------------------------------------\n";
  std::cout << "> orthoeop\n";
  std::cout << "----------------------------------------\n";

  const double result_ref[] = {-162.8386373279636530e0, 117.7907525842668974e0,
                               -23.39092370609808214e0};
  double result[3], result2[3];

  dso::datetime<dso::seconds> t(dso::modified_julian_day(47100), dso::seconds(0));
  iers2010::ortho_eop(t, result2[0], result2[1], result2[2]);

  iers2010::ortho_eop(47100e0, result[0], result[1], result[2]);
  for (int i = 0; i < 3; ++i) {
    printf("\nargs[%1d] = %12.6e microarcseconds", i,
           std::abs(result[i] - result_ref[i]));
    assert(approxEqual(result[i], result_ref[i]));
    assert(result[i] == result2[i]);
  }
  return 0;
}
