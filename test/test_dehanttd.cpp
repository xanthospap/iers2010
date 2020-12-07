#include "ggdatetime/dtcalendar.hpp"
#include "iers2010.hpp"
#include <iostream>

double input1[] = {
    4075578.385e0,       931852.890e0,        4801570.154e0,
    137859926952.015e0,  54228127881.4350e0,  23509422341.6960e0,
    -179996231.920342e0, -312468450.131567e0, -169288918.592160e0};
ngpt::datetime<ngpt::seconds> d1(ngpt::year(2009), ngpt::month(4),
                                 ngpt::day_of_month(13), ngpt::seconds(0));

double result1[] = {0.7700420357108125891e-01, 0.6304056321824967613e-01,
                    0.5516568152597246810e-01};

double input2[] = {1112189.660e0,       -4842955.026e0,     3985352.284e0,
                   -54537460436.2357e0, 130244288385.279e0, 56463429031.5996e0,
                   300396716.912e0,     243238281.451e0,    120548075.939e0};
ngpt::datetime<ngpt::seconds> d2(ngpt::year(2012), ngpt::month(7),
                                 ngpt::day_of_month(13), ngpt::seconds(0));

double result2[] = {-0.2036831479592075833e-01, 0.5658254776225972449e-01,
                    -0.7597679676871742227e-01};

double input3[] = {1112200.5696e0,      -4842957.8511e0,     3985345.9122e0,
                   100210282451.6279e0, 103055630398.3160e0, 56855096480.4475e0,
                   369817604.4348e0,    1897917.5258e0,      120804980.8284e0};
ngpt::datetime<ngpt::seconds> d3(ngpt::year(2015), ngpt::month(7),
                                 ngpt::day_of_month(15), ngpt::seconds(0));

double result3[] = {0.00509570869172363845e0, 0.0828663025983528700e0,
                    -0.0636634925404189617e0};

double input4[] = {1112152.8166e0,        3985496.1783e0,
                   8382471154.1312895e0,  10512408445.356153e0,
                   -5360583240.3763866e0, 380934092.93550891e0,
                   2871428.1904491195e0,  79015680.553570181e0};
ngpt::datetime<ngpt::seconds> d4(ngpt::year(2017), ngpt::month(1),
                                 ngpt::day_of_month(15), ngpt::seconds(0));

double result4[] = {0.0050957086917236384e0, 0.082866302598352870e0,
                    -0.063663492540418962e0};

int main() {
  double diff;
  double res[3];

  std::cout << "----------------------------------------\n";
  std::cout << "> dehanttideinel\n";
  std::cout << "----------------------------------------\n";

  std::cout << "/* Test a */\n";
  iers2010::dehanttideinel(input1, input1 + 3, input1 + 6, d1, res);
  for (int i = 0; i < 3; ++i) {
    diff = std::abs(res[i] - result1[i]);
    std::cout << "\tdiff = " << diff << "\n";
  }

  std::cout << "/* Test b */\n";
  iers2010::dehanttideinel(input2, input2 + 3, input2 + 6, d2, res);
  for (int i = 0; i < 3; ++i) {
    diff = std::abs(res[i] - result2[i]);
    std::cout << "\tdiff = " << diff << "\n";
  }

  std::cout << "/* Test c */\n";
  iers2010::dehanttideinel(input3, input3 + 3, input3 + 6, d3, res);
  for (int i = 0; i < 3; ++i) {
    diff = std::abs(res[i] - result3[i]);
    std::cout << "\tdiff = " << diff << "\n";
  }

  std::cout << "/* Test d */\n";
  iers2010::dehanttideinel(input4, input4 + 3, input4 + 6, d4, res);
  for (int i = 0; i < 3; ++i) {
    diff = std::abs(res[i] - result4[i]);
    std::cout << "\tdiff = " << diff << "\n";
  }

  return 0;
}
