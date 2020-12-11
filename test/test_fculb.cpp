#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

int main() {

  // testing fcul_b
  std::cout << "----------------------------------------\n";
  std::cout << "> fcul_b\n";
  std::cout << "----------------------------------------\n";

  const double result_ref = 3.800758725284345996e0;
  double f = iers2010::fcul_b(30.67166667e0, 2075e0, 224e0, 15e0);
#ifdef STRICT_TEST
  assert(approxEqual(f, result_ref));
#else
  printf("\nargs[%1d] = %12.6e ", 1, std::abs(f - result_ref));
  assert(std::abs(f - result_ref) < 1e-11);
#endif

  return 0;
}
