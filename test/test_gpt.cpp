#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

int main() {

  // testing gpt
  std::cout << "----------------------------------------\n";
  std::cout << "> gpt\n";
  std::cout << "----------------------------------------\n";

  const double result_ref[] = {918.0710638757363995e0, 19.31914181012882992e0,
                               -42.19185643717770517e0};
  const char *units[3] = {"hPa", "degrees Celsius", "meters"};
  double result[3];

  iers2010::gpt(55055e0, 0.6708665767e0, -1.393397187e0, 812.546e0, result[0],
                result[1], result[2]);
  for (int i = 0; i < 3; ++i) {
#ifdef STRICT_TEST
    assert(approxEqual(result[i], result_ref[i]));
#else
    printf("\nargs[%1d] = %12.6e %s", i, std::abs(result[i] - result_ref[i]),
           units[i]);
    // assert(std::abs(result[i] - result_ref[i]) < 1e-11);
#endif
  }
  return 0;
}
