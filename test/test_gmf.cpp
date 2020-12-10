#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

int main() {

  // testing gmf
  std::cout << "----------------------------------------\n";
  std::cout << "> gmf\n";
  std::cout << "----------------------------------------\n";

  const double result_ref[] = {3.425245519339138678e0, 3.449589116182419257e0};
  double result[2];

  iers2010::gmf(55055e0, 0.6708665767e0, -1.393397187e0, 844.715e0, 1.278564131e0, result[0], result[1]);
  for (int i = 0; i < 2; ++i) {
#ifdef STRICT_TEST
    assert(approxEqual(result[i], result_ref[i]));
#else
    printf("\nargs[%1d] = %12.6e radians", i,
           std::abs(result[i] - result_ref[i]));
    assert(std::abs(result[i] - result_ref[i]) < 1e-11);
#endif
  }
  return 0;
}
