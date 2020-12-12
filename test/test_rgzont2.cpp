#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

constexpr double _alg_accuracy_ = 1e-19;

int main() {

  // testing rg_zont2
  std::cout << "----------------------------------------\n";
  std::cout << "> rg_zont2\n";
  std::cout << "----------------------------------------\n";

  const char *units[] = {"seconds", "seconds / day", "radians / day"};
  const double result_ref[] = {7.983287678576557467E-002,
                               5.035331113978199288E-005,
                               -4.249711616463017E-014};
  double result[3];

  iers2010::rg_zont2(.07995893223819302e0, result[0], result[1], result[2]);
  for (int i = 0; i < 3; ++i) {
#ifdef STRICT_TEST
    assert(approxEqual(result[i], result_ref[i]));
#else
    printf("\nargs[%1d] = %12.6e %s", i, std::abs(result[i] - result_ref[i]),
           units[i]);
    assert(std::abs(result[i] - result_ref[i]) < _alg_accuracy_);
#endif
  }

  return 0;
}
