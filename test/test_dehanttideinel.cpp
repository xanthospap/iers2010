#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "iers2010.hpp"

//
// Test program for the function dehanttideinel
// To produce a new list of test cases, run the program makecc_test_cases.out
// found in fortran_impl/dehanttideinel (see readme there)
//

using iers2010::dehanttideinel;
using dso::Vector3;
using dso::datetime;
using dso::seconds;
using dso::year;
using dso::month;
using dso::day_of_month;

#define PRECISION 1e-9

TEST_CASE(
    "DEHANTTIDEINEL: Comparisson of results based on Fortran implementation.") {
  Vector3 xsta, xsun, xmon, xtide, res;
  datetime<seconds> t;

  SUBCASE("Example   1: Checking for descripancies  > PRECISION") {
    xsta.x() = 4075578.3850e0;
    xsta.y() = 931852.8900e0;
    xsta.z() = 4801570.1540e0;
    xsun.x() = 137859926952.0150e0;
    xsun.y() = 54228127881.4350e0;
    xsun.z() = 23509422341.6960e0;
    xmon.x() = -179996231.9203420e0;
    xmon.y() = -312468450.1315670e0;
    xmon.z() = -169288918.5921600e0;
    t = datetime<seconds>(year(2009), month(4), day_of_month(13), seconds(0.));
    xtide.x() = 0.0770042036e0;
    xtide.y() = 0.0630405632e0;
    xtide.z() = 0.0551656815e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example   2: Checking for descripancies  > PRECISION") {
    xsta.x() = 4075578.3850e0;
    xsta.y() = 931852.8900e0;
    xsta.z() = 4801570.1540e0;
    xsun.x() = 137859926952.0150e0;
    xsun.y() = 54228127881.4350e0;
    xsun.z() = 23509422341.6960e0;
    xmon.x() = -179996231.9203420e0;
    xmon.y() = -312468450.1315670e0;
    xmon.z() = -169288918.5921600e0;
    t = datetime<seconds>(year(2009), month(4), day_of_month(13),
                          seconds(3600.));
    xtide.x() = 0.0783698278e0;
    xtide.y() = 0.0634583620e0;
    xtide.z() = 0.0568293183e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example   3: Checking for descripancies  > PRECISION") {
    xsta.x() = 4075578.3850e0;
    xsta.y() = 931852.8900e0;
    xsta.z() = 4801570.1540e0;
    xsun.x() = 137859926952.0150e0;
    xsun.y() = 54228127881.4350e0;
    xsun.z() = 23509422341.6960e0;
    xmon.x() = -179996231.9203420e0;
    xmon.y() = -312468450.1315670e0;
    xmon.z() = -169288918.5921600e0;
    t = datetime<seconds>(year(2009), month(4), day_of_month(13),
                          seconds(14400.));
    xtide.x() = 0.0795930701e0;
    xtide.y() = 0.0641240473e0;
    xtide.z() = 0.0584000950e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example   4: Checking for descripancies  > PRECISION") {
    xsta.x() = 4075578.3850e0;
    xsta.y() = 931852.8900e0;
    xsta.z() = 4801570.1540e0;
    xsun.x() = 137859926952.0150e0;
    xsun.y() = 54228127881.4350e0;
    xsun.z() = 23509422341.6960e0;
    xmon.x() = -179996231.9203420e0;
    xmon.y() = -312468450.1315670e0;
    xmon.z() = -169288918.5921600e0;
    t = datetime<seconds>(year(2009), month(4), day_of_month(13),
                          seconds(25200.));
    xtide.x() = 0.0762814201e0;
    xtide.y() = 0.0636731386e0;
    xtide.z() = 0.0545293033e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example   5: Checking for descripancies  > PRECISION") {
    xsta.x() = 4075578.3850e0;
    xsta.y() = 931852.8900e0;
    xsta.z() = 4801570.1540e0;
    xsun.x() = 137859926952.0150e0;
    xsun.y() = 54228127881.4350e0;
    xsun.z() = 23509422341.6960e0;
    xmon.x() = -179996231.9203420e0;
    xmon.y() = -312468450.1315670e0;
    xmon.z() = -169288918.5921600e0;
    t = datetime<seconds>(year(2009), month(4), day_of_month(13),
                          seconds(36000.));
    xtide.x() = 0.0703658780e0;
    xtide.y() = 0.0623653692e0;
    xtide.z() = 0.0474836866e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example   6: Checking for descripancies  > PRECISION") {
    xsta.x() = 4075578.3850e0;
    xsta.y() = 931852.8900e0;
    xsta.z() = 4801570.1540e0;
    xsun.x() = 137859926952.0150e0;
    xsun.y() = 54228127881.4350e0;
    xsun.z() = 23509422341.6960e0;
    xmon.x() = -179996231.9203420e0;
    xmon.y() = -312468450.1315670e0;
    xmon.z() = -169288918.5921600e0;
    t = datetime<seconds>(year(2009), month(4), day_of_month(13),
                          seconds(46800.));
    xtide.x() = 0.0653134822e0;
    xtide.y() = 0.0609658655e0;
    xtide.z() = 0.0414028296e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example   7: Checking for descripancies  > PRECISION") {
    xsta.x() = 4075578.3850e0;
    xsta.y() = 931852.8900e0;
    xsta.z() = 4801570.1540e0;
    xsun.x() = 137859926952.0150e0;
    xsun.y() = 54228127881.4350e0;
    xsun.z() = 23509422341.6960e0;
    xmon.x() = -179996231.9203420e0;
    xmon.y() = -312468450.1315670e0;
    xmon.z() = -169288918.5921600e0;
    t = datetime<seconds>(year(2009), month(4), day_of_month(13),
                          seconds(57600.));
    xtide.x() = 0.0640927782e0;
    xtide.y() = 0.0602968372e0;
    xtide.z() = 0.0398698701e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example   8: Checking for descripancies  > PRECISION") {
    xsta.x() = 4075578.3850e0;
    xsta.y() = 931852.8900e0;
    xsta.z() = 4801570.1540e0;
    xsun.x() = 137859926952.0150e0;
    xsun.y() = 54228127881.4350e0;
    xsun.z() = 23509422341.6960e0;
    xmon.x() = -179996231.9203420e0;
    xmon.y() = -312468450.1315670e0;
    xmon.z() = -169288918.5921600e0;
    t = datetime<seconds>(year(2009), month(4), day_of_month(13),
                          seconds(68400.));
    xtide.x() = 0.0674276382e0;
    xtide.y() = 0.0607540175e0;
    xtide.z() = 0.0438035685e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example   9: Checking for descripancies  > PRECISION") {
    xsta.x() = 4075578.3850e0;
    xsta.y() = 931852.8900e0;
    xsta.z() = 4801570.1540e0;
    xsun.x() = 137859926952.0150e0;
    xsun.y() = 54228127881.4350e0;
    xsun.z() = 23509422341.6960e0;
    xmon.x() = -179996231.9203420e0;
    xmon.y() = -312468450.1315670e0;
    xmon.z() = -169288918.5921600e0;
    t = datetime<seconds>(year(2009), month(4), day_of_month(13),
                          seconds(79200.));
    xtide.x() = 0.0733669187e0;
    xtide.y() = 0.0620723474e0;
    xtide.z() = 0.0509121860e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  10: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112189.6600e0;
    xsta.y() = -4842955.0260e0;
    xsta.z() = 3985352.2840e0;
    xsun.x() = -54537460436.2357e0;
    xsun.y() = 130244288385.2790e0;
    xsun.z() = 56463429031.5996e0;
    xmon.x() = 300396716.9120000e0;
    xmon.y() = 243238281.4510000e0;
    xmon.z() = 120548075.9390000e0;
    t = datetime<seconds>(year(2012), month(7), day_of_month(13), seconds(0.));
    xtide.x() = -0.0203683148e0;
    xtide.y() = 0.0565825478e0;
    xtide.z() = -0.0759767968e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  11: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112189.6600e0;
    xsta.y() = -4842955.0260e0;
    xsta.z() = 3985352.2840e0;
    xsun.x() = -54537460436.2357e0;
    xsun.y() = 130244288385.2790e0;
    xsun.z() = 56463429031.5996e0;
    xmon.x() = 300396716.9120000e0;
    xmon.y() = 243238281.4510000e0;
    xmon.z() = 120548075.9390000e0;
    t = datetime<seconds>(year(2012), month(7), day_of_month(13),
                          seconds(7200.));
    xtide.x() = -0.0194853460e0;
    xtide.y() = 0.0535951404e0;
    xtide.z() = -0.0735316514e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  12: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112189.6600e0;
    xsta.y() = -4842955.0260e0;
    xsta.z() = 3985352.2840e0;
    xsun.x() = -54537460436.2357e0;
    xsun.y() = 130244288385.2790e0;
    xsun.z() = 56463429031.5996e0;
    xmon.x() = 300396716.9120000e0;
    xmon.y() = 243238281.4510000e0;
    xmon.z() = 120548075.9390000e0;
    t = datetime<seconds>(year(2012), month(7), day_of_month(13),
                          seconds(18000.));
    xtide.x() = -0.0191768616e0;
    xtide.y() = 0.0537543116e0;
    xtide.z() = -0.0735783058e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  13: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112189.6600e0;
    xsta.y() = -4842955.0260e0;
    xsta.z() = 3985352.2840e0;
    xsun.x() = -54537460436.2357e0;
    xsun.y() = 130244288385.2790e0;
    xsun.z() = 56463429031.5996e0;
    xmon.x() = 300396716.9120000e0;
    xmon.y() = 243238281.4510000e0;
    xmon.z() = 120548075.9390000e0;
    t = datetime<seconds>(year(2012), month(7), day_of_month(13),
                          seconds(28800.));
    xtide.x() = -0.0201951653e0;
    xtide.y() = 0.0591565198e0;
    xtide.z() = -0.0778592381e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  14: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112189.6600e0;
    xsta.y() = -4842955.0260e0;
    xsta.z() = 3985352.2840e0;
    xsun.x() = -54537460436.2357e0;
    xsun.y() = 130244288385.2790e0;
    xsun.z() = 56463429031.5996e0;
    xmon.x() = 300396716.9120000e0;
    xmon.y() = 243238281.4510000e0;
    xmon.z() = 120548075.9390000e0;
    t = datetime<seconds>(year(2012), month(7), day_of_month(13),
                          seconds(39600.));
    xtide.x() = -0.0219442018e0;
    xtide.y() = 0.0666352730e0;
    xtide.z() = -0.0838602295e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  15: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112189.6600e0;
    xsta.y() = -4842955.0260e0;
    xsta.z() = 3985352.2840e0;
    xsun.x() = -54537460436.2357e0;
    xsun.y() = 130244288385.2790e0;
    xsun.z() = 56463429031.5996e0;
    xmon.x() = 300396716.9120000e0;
    xmon.y() = 243238281.4510000e0;
    xmon.z() = 120548075.9390000e0;
    t = datetime<seconds>(year(2012), month(7), day_of_month(13),
                          seconds(50400.));
    xtide.x() = -0.0233954682e0;
    xtide.y() = 0.0717939353e0;
    xtide.z() = -0.0880480580e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  16: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112189.6600e0;
    xsta.y() = -4842955.0260e0;
    xsta.z() = 3985352.2840e0;
    xsun.x() = -54537460436.2357e0;
    xsun.y() = 130244288385.2790e0;
    xsun.z() = 56463429031.5996e0;
    xmon.x() = 300396716.9120000e0;
    xmon.y() = 243238281.4510000e0;
    xmon.z() = 120548075.9390000e0;
    t = datetime<seconds>(year(2012), month(7), day_of_month(13),
                          seconds(61200.));
    xtide.x() = -0.0236931928e0;
    xtide.y() = 0.0715920306e0;
    xtide.z() = -0.0879487649e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  17: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112189.6600e0;
    xsta.y() = -4842955.0260e0;
    xsta.z() = 3985352.2840e0;
    xsun.x() = -54537460436.2357e0;
    xsun.y() = 130244288385.2790e0;
    xsun.z() = 56463429031.5996e0;
    xmon.x() = 300396716.9120000e0;
    xmon.y() = 243238281.4510000e0;
    xmon.z() = 120548075.9390000e0;
    t = datetime<seconds>(year(2012), month(7), day_of_month(13),
                          seconds(72000.));
    xtide.x() = -0.0226591193e0;
    xtide.y() = 0.0661380308e0;
    xtide.z() = -0.0836064481e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  18: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112189.6600e0;
    xsta.y() = -4842955.0260e0;
    xsta.z() = 3985352.2840e0;
    xsun.x() = -54537460436.2357e0;
    xsun.y() = 130244288385.2790e0;
    xsun.z() = 56463429031.5996e0;
    xmon.x() = 300396716.9120000e0;
    xmon.y() = 243238281.4510000e0;
    xmon.z() = 120548075.9390000e0;
    t = datetime<seconds>(year(2012), month(7), day_of_month(13),
                          seconds(82800.));
    xtide.x() = -0.0208992193e0;
    xtide.y() = 0.0586318012e0;
    xtide.z() = -0.0775626218e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  19: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15), seconds(0.));
    xtide.x() = 0.0050957087e0;
    xtide.y() = 0.0828663026e0;
    xtide.z() = -0.0636634925e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  20: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15), seconds(0.));
    xtide.x() = 0.0050957087e0;
    xtide.y() = 0.0828663026e0;
    xtide.z() = -0.0636634925e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  21: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(150.));
    xtide.x() = 0.0051186773e0;
    xtide.y() = 0.0827794975e0;
    xtide.z() = -0.0635929941e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  22: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(300.));
    xtide.x() = 0.0051415541e0;
    xtide.y() = 0.0826932610e0;
    xtide.z() = -0.0635229422e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  23: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(450.));
    xtide.x() = 0.0051643365e0;
    xtide.y() = 0.0826076034e0;
    xtide.z() = -0.0634533452e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  24: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(600.));
    xtide.x() = 0.0051870217e0;
    xtide.y() = 0.0825225349e0;
    xtide.z() = -0.0633842115e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  25: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(750.));
    xtide.x() = 0.0052096071e0;
    xtide.y() = 0.0824380655e0;
    xtide.z() = -0.0633155491e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  26: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(900.));
    xtide.x() = 0.0052320899e0;
    xtide.y() = 0.0823542055e0;
    xtide.z() = -0.0632473664e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  27: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(1050.));
    xtide.x() = 0.0052544674e0;
    xtide.y() = 0.0822709646e0;
    xtide.z() = -0.0631796714e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  28: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(1200.));
    xtide.x() = 0.0052767370e0;
    xtide.y() = 0.0821883529e0;
    xtide.z() = -0.0631124721e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  29: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(1350.));
    xtide.x() = 0.0052988961e0;
    xtide.y() = 0.0821063802e0;
    xtide.z() = -0.0630457765e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  30: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(1500.));
    xtide.x() = 0.0053209420e0;
    xtide.y() = 0.0820250562e0;
    xtide.z() = -0.0629795925e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  31: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(1650.));
    xtide.x() = 0.0053428721e0;
    xtide.y() = 0.0819443906e0;
    xtide.z() = -0.0629139281e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  32: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(1800.));
    xtide.x() = 0.0053646838e0;
    xtide.y() = 0.0818643931e0;
    xtide.z() = -0.0628487911e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  33: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(1950.));
    xtide.x() = 0.0053863745e0;
    xtide.y() = 0.0817850730e0;
    xtide.z() = -0.0627841890e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  34: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(2100.));
    xtide.x() = 0.0054079416e0;
    xtide.y() = 0.0817064400e0;
    xtide.z() = -0.0627201298e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  35: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(2250.));
    xtide.x() = 0.0054293826e0;
    xtide.y() = 0.0816285032e0;
    xtide.z() = -0.0626566209e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  36: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(2400.));
    xtide.x() = 0.0054506948e0;
    xtide.y() = 0.0815512721e0;
    xtide.z() = -0.0625936699e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  37: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(2550.));
    xtide.x() = 0.0054718758e0;
    xtide.y() = 0.0814747557e0;
    xtide.z() = -0.0625312843e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  38: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(2700.));
    xtide.x() = 0.0054929231e0;
    xtide.y() = 0.0813989633e0;
    xtide.z() = -0.0624694715e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  39: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(2850.));
    xtide.x() = 0.0055138340e0;
    xtide.y() = 0.0813239037e0;
    xtide.z() = -0.0624082388e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  40: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(3000.));
    xtide.x() = 0.0055346063e0;
    xtide.y() = 0.0812495860e0;
    xtide.z() = -0.0623475936e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  41: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(3150.));
    xtide.x() = 0.0055552373e0;
    xtide.y() = 0.0811760190e0;
    xtide.z() = -0.0622875430e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  42: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(3300.));
    xtide.x() = 0.0055757246e0;
    xtide.y() = 0.0811032114e0;
    xtide.z() = -0.0622280941e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  43: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112200.5696e0;
    xsta.y() = -4842957.8511e0;
    xsta.z() = 3985345.9122e0;
    xsun.x() = 100210282451.6279e0;
    xsun.y() = 103055630398.3160e0;
    xsun.z() = 56855096480.4475e0;
    xmon.x() = 369817604.4348000e0;
    xmon.y() = 1897917.5258000e0;
    xmon.z() = 120804980.8284000e0;
    t = datetime<seconds>(year(2015), month(7), day_of_month(15),
                          seconds(3450.));
    xtide.x() = 0.0055960659e0;
    xtide.y() = 0.0810311719e0;
    xtide.z() = -0.0621692541e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  44: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15), seconds(0.));
    xtide.x() = -18.2173575819e0;
    xtide.y() = -23.5053483765e0;
    xtide.z() = 12.0976113822e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  45: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15), seconds(0.));
    xtide.x() = -18.2173575819e0;
    xtide.y() = -23.5053483765e0;
    xtide.z() = 12.0976113822e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  46: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(150.));
    xtide.x() = -18.2173796113e0;
    xtide.y() = -23.5052634891e0;
    xtide.z() = 12.0975426898e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  47: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(300.));
    xtide.x() = -18.2174015874e0;
    xtide.y() = -23.5051790139e0;
    xtide.z() = 12.0974743173e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  48: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(450.));
    xtide.x() = -18.2174235076e0;
    xtide.y() = -23.5050949611e0;
    xtide.z() = 12.0974062726e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  49: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(600.));
    xtide.x() = -18.2174453692e0;
    xtide.y() = -23.5050113406e0;
    xtide.z() = 12.0973385640e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  50: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(750.));
    xtide.x() = -18.2174671697e0;
    xtide.y() = -23.5049281625e0;
    xtide.z() = 12.0972711995e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  51: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(900.));
    xtide.x() = -18.2174889064e0;
    xtide.y() = -23.5048454366e0;
    xtide.z() = 12.0972041870e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  52: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(1050.));
    xtide.x() = -18.2175105768e0;
    xtide.y() = -23.5047631728e0;
    xtide.z() = 12.0971375347e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  53: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(1200.));
    xtide.x() = -18.2175321783e0;
    xtide.y() = -23.5046813809e0;
    xtide.z() = 12.0970712505e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  54: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(1350.));
    xtide.x() = -18.2175537083e0;
    xtide.y() = -23.5046000706e0;
    xtide.z() = 12.0970053421e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  55: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(1500.));
    xtide.x() = -18.2175751643e0;
    xtide.y() = -23.5045192517e0;
    xtide.z() = 12.0969398176e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  56: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(1650.));
    xtide.x() = -18.2175965436e0;
    xtide.y() = -23.5044389338e0;
    xtide.z() = 12.0968746847e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  57: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(1800.));
    xtide.x() = -18.2176178438e0;
    xtide.y() = -23.5043591264e0;
    xtide.z() = 12.0968099511e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  58: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(1950.));
    xtide.x() = -18.2176390623e0;
    xtide.y() = -23.5042798391e0;
    xtide.z() = 12.0967456246e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  59: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(2100.));
    xtide.x() = -18.2176601966e0;
    xtide.y() = -23.5042010813e0;
    xtide.z() = 12.0966817128e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  60: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(2250.));
    xtide.x() = -18.2176812442e0;
    xtide.y() = -23.5041228623e0;
    xtide.z() = 12.0966182233e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  61: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(2400.));
    xtide.x() = -18.2177022025e0;
    xtide.y() = -23.5040451915e0;
    xtide.z() = 12.0965551638e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  62: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(2550.));
    xtide.x() = -18.2177230692e0;
    xtide.y() = -23.5039680782e0;
    xtide.z() = 12.0964925416e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  63: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(2700.));
    xtide.x() = -18.2177438416e0;
    xtide.y() = -23.5038915315e0;
    xtide.z() = 12.0964303643e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  64: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(2850.));
    xtide.x() = -18.2177645173e0;
    xtide.y() = -23.5038155605e0;
    xtide.z() = 12.0963686392e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  65: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(3000.));
    xtide.x() = -18.2177850939e0;
    xtide.y() = -23.5037401743e0;
    xtide.z() = 12.0963073737e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  66: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(3150.));
    xtide.x() = -18.2178055689e0;
    xtide.y() = -23.5036653819e0;
    xtide.z() = 12.0962465751e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  67: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(3300.));
    xtide.x() = -18.2178259399e0;
    xtide.y() = -23.5035911921e0;
    xtide.z() = 12.0961862506e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
  SUBCASE("Example  68: Checking for descripancies  > PRECISION") {
    xsta.x() = 1112152.8166e0;
    xsta.y() = -4842857.5435e0;
    xsta.z() = 3985496.1783e0;
    xsun.x() = 8382471154.1313e0;
    xsun.y() = 10512408445.3562e0;
    xsun.z() = -5360583240.3764e0;
    xmon.x() = 380934092.9355089e0;
    xmon.y() = 2871428.1904491e0;
    xmon.z() = 79015680.5535702e0;
    t = datetime<seconds>(year(2017), month(1), day_of_month(15),
                          seconds(3450.));
    xtide.x() = -18.2178462044e0;
    xtide.y() = -23.5035176138e0;
    xtide.z() = 12.0961264075e0;
    res = dehanttideinel(xsta, xsun, xmon, t);
    REQUIRE(std::abs(res.x() - xtide.x()) < PRECISION);
    REQUIRE(std::abs(res.y() - xtide.y()) < PRECISION);
    REQUIRE(std::abs(res.z() - xtide.z()) < PRECISION);
  }
}