#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

int main() {

  // testing fundarg
  std::cout << "----------------------------------------\n";
  std::cout << "> fcnnut\n";
  std::cout << "----------------------------------------\n";

  const double result_ref[] = {
    -176.8012290066270680e0,
    -93.51855308903756736e0,
    3.745573770491803067e0,
    3.745573770491803067e0};
  double result[4];

  iers2010::fcnnut(54790e0, result[0], result[1], result[2], result[3]);
  for (int i = 0; i < 4; ++i) {
#ifdef STRICT_TEST
    assert(approxEqual(result[i], result_ref[i]));
#else
    printf("\nargs[%1d] = %12.6e microarcseconds", i, std::abs(result[i]- result_ref[i]));
    assert(std::abs(result[i]- result_ref[i])<1e-11);
#endif
  }
  return 0;
}
