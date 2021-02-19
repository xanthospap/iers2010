#include "ggdatetime/dtfund.hpp"
#include "hardisp.hpp"
#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

constexpr double _alg_accuracy_ = 1e-8;

int main() {

  // testing tdfrph
  std::cout << "----------------------------------------\n";
  std::cout << "> tdfrph\n";
  std::cout << "----------------------------------------\n";

  auto t1 =
      ngpt::datetime<ngpt::seconds>{ngpt::year(2009), ngpt::day_of_year(176)};
  auto t2 = ngpt::datetime<ngpt::seconds>{
      ngpt::year(2009), ngpt::day_of_year(176), ngpt::hours(12),
      ngpt::minutes(1), ngpt::seconds(45)};
  const ngpt::datetime<ngpt::seconds> epochs[] = {t1, t2};
  const int idood[2][6] = {{2, 0, 0, 0, 0, 0}, {2, 0, 0, 0, 0, 0}};
  const double result_ref[2][2] = {{1.93227361605688e0, 303.343338720297e0},
                                   {1.93227361605689e0, 291.997959322689e0}};
  double result[2];

  for (int k = 0; k < 2; k++) {
    iers2010::hisp::tdfrph(idood[k], epochs[k], result[0], result[1]);
    for (int i = 0; i < 2; ++i) {
#ifdef STRICT_TEST
      assert(approxEqual(result[i], result_ref[k][i]));
#else
      printf("\nargs[%1d] = %12.6e ", i,
             std::abs(result[i] - result_ref[k][i]));
      assert(std::abs(result[i] - result_ref[k][i]) < _alg_accuracy_);
#endif
    }
  }
  return 0;
}
