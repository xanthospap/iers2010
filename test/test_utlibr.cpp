#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

constexpr double _alg_accuracy_ = 1e-6;

int main() {

  // testing utlibr
  std::cout << "----------------------------------------\n";
  std::cout << "> utlibr\n";
  std::cout << "----------------------------------------\n";

  const double results[] = {2.441143834386761746e0, // test case A
                            -14.78971247349449492e0,
                            -2.655705844335680244e0, // test case B
                            27.39445826599846967e0};
  double dut1, dlod;

  iers2010::utlibr(44239.1e0, dut1, dlod);
#ifdef STRICT_TEST
  assert(approxEqual(dut1, results[0]));
  assert(approxEqual(dlod, results[1]));
#else
  printf("\ndut1= %12.6e mas", std::abs(dut1 - results[0]));
  printf("\ndlod= %12.6e mas / day", std::abs(dlod - results[1]));
  assert(std::abs(dut1 - results[0]) < _alg_accuracy_ &&
         std::abs(dlod - results[1]) < _alg_accuracy_);
#endif

  dso::datetime<dso::seconds> t(dso::year(1980), dso::month(1),
                                dso::day_of_month(1), dso::hours(2),
                                dso::minutes(24), dso::seconds(0));
  double dut1_, dlod_;
  iers2010::utlibr(t, dut1_, dlod_);
  assert(dut1 == dut1_ && dlod == dlod_);

  iers2010::utlibr(55227.4e0, dut1, dlod);
#ifdef STRICT_TEST
  assert(approxEqual(dut1, results[2]));
  assert(approxEqual(dlod, results[3]));
#else
  printf("\ndut1= %12.6e mas", std::abs(dut1 - results[2]));
  printf("\ndlod= %12.6e mas / day", std::abs(dlod - results[3]));
  assert(std::abs(dut1 - results[2]) < _alg_accuracy_ &&
         std::abs(dlod - results[3]) < _alg_accuracy_);
#endif
  
  return 0;
}
