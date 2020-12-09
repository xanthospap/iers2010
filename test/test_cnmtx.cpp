#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

int main() {

  // testing cnmtx
  std::cout << "----------------------------------------\n";
  std::cout << "> cnmtx\n";
  std::cout << "----------------------------------------\n";

  const double result_ref[] = {
15.35873641938967360e0,
9.784941251812741214e0,
-5.520740128266865554e0,
3.575314211234633888e0,
-13.93717453496387648e0,
-9.167400321705855504e0,
5.532815475865292321e0,
9.558741883500834646e0,
-10.22541212627272600e0,
0.8367570529461261231e0,
1.946355176475630611e0,
-13.55702062247304696e0};
  double result[12];

  iers2010::oeop::cnmtx(54964.0e0, result);
  for (int i = 0; i < 12; ++i) {
#ifdef STRICT_TEST
    assert(approxEqual(result[i], result_ref[i]));
#else
    printf("\nargs[%1d] = %12.6e ", i,
           std::abs(result[i] - result_ref[i]));
    assert(std::abs(result[i] - result_ref[i]) < 1e-11);
#endif
  }
  return 0;
}
