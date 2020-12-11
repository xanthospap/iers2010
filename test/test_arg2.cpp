#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

int main() {

  // testing arg2
  std::cout << "----------------------------------------\n";
  std::cout << "> arg2\n";
  std::cout << "----------------------------------------\n";

  const double result_ref[] = {
      2.849663065753787805e0, 6.28318080000000023e0,  4.926040134021299366e0,
      1.608450491115348768e0, 2.375021572352622456e0, 0.4746414933980958040e0,
      3.908159227647345801e0, 2.551018561669245344e0, 5.041990012540757959e0,
      4.206816878908014701e0, 1.608463638294885811e0,
  };
  double result[11];

  iers2010::arg2(2008, 311.5e0, result);
  for (int i = 0; i < 11; ++i) {
#ifdef STRICT_TEST
    assert(approxEqual(result[i], result_ref[i]));
#else
    printf("\nargs[%1d] = %12.6e radians", i,
           std::abs(result[i] - result_ref[i]));
    // assert(std::abs(result[i] - result_ref[i]) < 1e-11);
#endif
  }
  return 0;
}
