#include "ggdatetime/dtcalendar.hpp"
#include "test_help.hpp"
#include <cassert>
#include "iers2010.hpp"
#include <iostream>
#include <map>

const double input1[] = {
    4075578.385e0,       931852.890e0,        4801570.154e0,
    137859926952.015e0,  54228127881.4350e0,  23509422341.6960e0,
    -179996231.920342e0, -312468450.131567e0, -169288918.592160e0};
const ngpt::datetime<ngpt::seconds> d1(ngpt::year(2009), ngpt::month(4),
                                 ngpt::day_of_month(13), ngpt::seconds(0));

const double result1[] = {0.7700420357108125891e-01, 0.6304056321824967613e-01,
                    0.5516568152597246810e-01};

const double input2[] = {1112189.660e0,       -4842955.026e0,     3985352.284e0,
                   -54537460436.2357e0, 130244288385.279e0, 56463429031.5996e0,
                   300396716.912e0,     243238281.451e0,    120548075.939e0};
const ngpt::datetime<ngpt::seconds> d2(ngpt::year(2012), ngpt::month(7),
                                 ngpt::day_of_month(13), ngpt::seconds(0));

const double result2[] = {-0.2036831479592075833e-01, 0.5658254776225972449e-01,
                    -0.7597679676871742227e-01};

const double input3[] = {1112200.5696e0,      -4842957.8511e0,     3985345.9122e0,
                   100210282451.6279e0, 103055630398.3160e0, 56855096480.4475e0,
                   369817604.4348e0,    1897917.5258e0,      120804980.8284e0};
const ngpt::datetime<ngpt::seconds> d3(ngpt::year(2015), ngpt::month(7),
                                 ngpt::day_of_month(15), ngpt::seconds(0));

const double result3[] = {0.00509570869172363845e0, 0.0828663025983528700e0,
                    -0.0636634925404189617e0};

const double input4[] = {1112152.8166e0,        3985496.1783e0,
                   8382471154.1312895e0,  10512408445.356153e0,
                   -5360583240.3763866e0, 380934092.93550891e0,
                   2871428.1904491195e0,  79015680.553570181e0};
const ngpt::datetime<ngpt::seconds> d4(ngpt::year(2017), ngpt::month(1),
                                 ngpt::day_of_month(15), ngpt::seconds(0));

const double result4[] = {0.0050957086917236384e0, 0.082866302598352870e0,
                    -0.063663492540418962e0};

struct testCaseData {
  const ngpt::datetime<ngpt::seconds> t_;
  const double *result_, *input_;
  const double* xsta() const { return input_; }
  const double* xsun() const { return input_+3; }
  const double* xmon() const { return input_+6; }
};

typedef std::map<char, testCaseData> test_case_map;
test_case_map tests = {{'a', {d1, result1, input1}}, 
{'b', {d2, result2, input2}}, 
{'c', {d3, result3, input3}}, 
{'d', {d4, result4, input4}}}; 


int main() {

  std::cout << "----------------------------------------\n";
  std::cout << "> dehanttideinel\n";
  std::cout << "----------------------------------------\n";
  
  double res[3];
  for (const auto& test : tests) {
    printf("\n// test case %c", test.first);
    iers2010::dehanttideinel(test.second.xsta(), test.second.xsun(), test.second.xmon(), test.second.t_, res);
    for (int i = 0; i < 3; ++i) {
#ifdef STRICT_TEST
      assert(approxEqual(res[i], test.second.input_[i]));
#else
      printf("\nargs[%1d] = %12.6e meters", i, std::abs(res[i]- test.second.result_[i]));
      //assert(std::abs(res[i]- test.second.input_[i])<1e-11);
#endif
    }
  }

  return 0;
}
