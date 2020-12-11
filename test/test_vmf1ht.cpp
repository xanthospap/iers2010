#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

int main() {

  // testing vmf1_ht
  std::cout << "----------------------------------------\n";
  std::cout << "> vmf1_ht\n";
  std::cout << "----------------------------------------\n";

  const double result_ref[] = {3.425088087972572470e0, 3.448299714692572238e0};
  double result[2];

  iers2010::vmf1_ht(0.00127683e0, 0.00060955e0, 55055e0, 0.6708665767e0,
                    824.17e0, 1.278564131e0, result[0], result[1]);
  for (int i = 0; i < 2; ++i) {
#ifdef STRICT_TEST
    assert(approxEqual(result[i], result_ref[i]));
#else
    printf("\nargs[%1d] = %12.6e ", i, std::abs(result[i] - result_ref[i]));
    // assert(std::abs(result[i] - result_ref[i]) < 1e-11);
#endif
  }
  return 0;
}
