#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

constexpr double _alg_accuray_ = 1e-5;

int main() {

  // testing fcul_zd_hpa
  std::cout << "----------------------------------------\n";
  std::cout << "> fcul_zd_hpa\n";
  std::cout << "----------------------------------------\n";

  const double result_ref[] = {1.935225924846803114e0, 1.932992176591644462e0,
                               0.2233748255158703871e-02};
  double result[3];

  iers2010::fcul_zd_hpa(30.67166667e0, 2010.344e0, 798.4188e0, 14.322e0,
                        0.532e0, result[0], result[1], result[2]);
  for (int i = 0; i < 3; i++) {
#ifdef STRICT_TEST
    assert(approxEqual(result[i], result_ref[i]));
#else
    printf("\nargs[%1d] = %12.6e meters", 1,
           std::abs(result[i] - result_ref[i]));
    assert(std::abs(result[i] - result_ref[i]) < _alg_accuray_);
#endif
  }
  return 0;
}
