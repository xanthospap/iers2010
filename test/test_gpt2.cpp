#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

#define pi 3.1415926535897932384626433e0

int main() {

  // testing gpt2
  std::cout << "----------------------------------------\n";
  std::cout << "> gpt2\n";
  std::cout << "----------------------------------------\n";

  const double result_ref[][7] = {
      {1002.56e0, 22.12e0, -6.53e0, 15.63e0, 0.0012647e0, 0.0005726e0, 44.06e0},
      {1003.49e0, 11.95e0, -5.47e0, 9.58e0, 0.0012395e0, 0.0005560e0, 44.06e0}};
  const char *units[] = {"hPa", "degrees Celsius", "deg/km", "hPa", "",
                         "",    "meters"};
  const double input[][6] = {
      {56141.e0, 48.20e0 * pi / 180.e0, 16.37e0 * pi / 180.e0, 156.e0, 1, 0},
      {56141.e0, 48.20e0 * pi / 180.e0, 16.37e0 * pi / 180.e0, 156.e0, 1, 1}};

  const double alg_accuracy_[] = {1e-2, 1e-2, 1e-2, 1e-2, 1e-7, 1e-7, 1e-2};
  double result[7];

  for (int t = 0; t < 2; t++) {
    iers2010::gpt2(input[t][0], &input[t][1], &input[t][2], &input[t][3],
                   (int)input[t][4], (int)input[t][5], &result[0], &result[1],
                   &result[2], &result[3], &result[4], &result[5], &result[6]);
    for (int i = 0; i < 7; ++i) {
#ifdef STRICT_TEST
      assert(approxEqual(result[i], result_ref[t][i]));
#else
      printf("\nargs[%1d] = %12.6e %s", i,
             std::abs(result[i] - result_ref[t][i]), units[i]);
      assert(std::abs(result[i] - result_ref[t][i]) < alg_accuracy_[i]);
#endif
    }
  }
  return 0;
}
