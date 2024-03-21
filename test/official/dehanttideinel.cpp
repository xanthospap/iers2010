#include "fundarg.hpp"
#include "solid_earth_tide.hpp"
#include <cstdio>

using namespace dso;

int main() {

  dso::CartesianCrd rsta, rmon, rsun;
  SolidEarthTide set(3.986004418e14, 6378136.6e0, 4.9048695e12,
                     1.32712442099e20);

  /* Example 1 */
  // printf("\nTest Case 1:\n");
  rsta.x() = 4075578.385e0;       //[m]
  rsta.y() = 931852.890e0;        //[m]
  rsta.z() = 4801570.154e0;       //[m]
  rsun.x() = 137859926952.015e0;  //[m]
  rsun.y() = 54228127881.4350e0;  //[m]
  rsun.z() = 23509422341.6960e0;  //[m]
  rmon.x() = -179996231.920342e0; //[m]
  rmon.y() = -312468450.131567e0; //[m]
  rmon.z() = -169288918.592160e0; //[m]

  TwoPartDateUTC utc(dso::datetime<nanoseconds>(
      year(2009), month(4), day_of_month(13), nanoseconds(0)));
  MjdEpoch tt(utc.utc2tt());
  MjdEpoch ut(tt.tt2ut1(0.3089055e0));

  /* Fundamental arguments (IERS 2003) */
  double fa[14] = {
      fal03(tt),  falp03(tt), faf03(tt),  fad03(tt),  faom03(tt),
      fame03(tt), fave03(tt), fae03(tt),  fama03(tt), faju03(tt),
      fasa03(tt), faur03(tt), fane03(tt), fapa03(tt),
  };

  auto dr = set.displacement(tt, ut, rsta.mv, rmon.mv, rsun.mv, fa);
  //printf("Displacement : %+.6f %+.6f %+.6f\n", dr(0), dr(1), dr(2));
  //printf("Diffs        : %+.4f %+.4f %+.4f\n",
  //       dr(0) - 0.7700420357108125891e-01, dr(1) - 0.6304056321824967613e-01,
  //       dr(2) - 0.5516568152597246810e-01);
  
  /* Example  2 */
  printf("\nTest Case 2:\n");
  rsta.x() = 1112189.660e0;//[m]
  rsta.y() = -4842955.026e0;//[m]
  rsta.z() = 3985352.284e0;//[m]
  rsun.x() = -54537460436.2357e0;//[m]
  rsun.y() = 130244288385.279e0;//[m]
  rsun.z() = 56463429031.5996e0;//[m]
  rmon.x() = 300396716.912e0;//[m]
  rmon.y() = 243238281.451e0;//[m]
  rmon.z() = 120548075.939e0;//[m]

  utc = TwoPartDateUTC(dso::datetime<nanoseconds>(
      year(2012), month(7), day_of_month(13), nanoseconds(0)));
  tt=(utc.utc2tt());
  ut=(tt.tt2ut1(0.4128669));

  /* Fundamental arguments (IERS 2003) */
  fa[0] = fal03(tt);
  fa[1] = falp03(tt);
  fa[2] = faf03(tt);
  fa[3] = fad03(tt);
  fa[4] = faom03(tt);
  fa[5] = fame03(tt);
  fa[6] = fave03(tt);
  fa[7] = fae03(tt);
  fa[8] = fama03(tt);
  fa[9] = faju03(tt);
  fa[10] = fasa03(tt);
  fa[11] = faur03(tt);
  fa[12] = fane03(tt);
  fa[13] = fapa03(tt);

  dr = set.displacement(tt, ut, rsta.mv, rmon.mv, rsun.mv, fa);
  //printf("Displacement : %+.6f %+.6f %+.6f\n", dr(0), dr(1), dr(2));
  //printf("Diffs        : %+.4f %+.4f %+.4f\n",
  //       dr(0) - (-0.2036831479592075833e-01),
  //       dr(1) - 0.5658254776225972449e-01,
  //       dr(2) - (-0.7597679676871742227e-01));
  
  /* Example  3 */
  printf("\nTest Case 3:\n");
  rsta.x() = 1112200.5696e0;
  rsta.y() =  -4842957.8511e0;
  rsta.z() = 3985345.9122e0;
  rsun.x() = 100210282451.6279e0;
  rsun.y() =103055630398.3160e0;
  rsun.z() =56855096480.4475e0;
  rmon.x() = 369817604.4348e0;
  rmon.y() = 1897917.5258e0;
  rmon.z() = 120804980.8284e0;

  utc = TwoPartDateUTC(dso::datetime<nanoseconds>(
      year(2015), month(7), day_of_month(15), nanoseconds(0)));
  tt=(utc.utc2tt());
  ut=(tt.tt2ut1(0.3101196));

  /* Fundamental arguments (IERS 2003) */
  fa[0] = fal03(tt);
  fa[1] = falp03(tt);
  fa[2] = faf03(tt);
  fa[3] = fad03(tt);
  fa[4] = faom03(tt);
  fa[5] = fame03(tt);
  fa[6] = fave03(tt);
  fa[7] = fae03(tt);
  fa[8] = fama03(tt);
  fa[9] = faju03(tt);
  fa[10] = fasa03(tt);
  fa[11] = faur03(tt);
  fa[12] = fane03(tt);
  fa[13] = fapa03(tt);

  dr = set.displacement(tt, ut, rsta.mv, rmon.mv, rsun.mv, fa);
  //printf("Displacement : %+.6f %+.6f %+.6f\n", dr(0), dr(1), dr(2));
  //printf("Diffs        : %+.4f %+.4f %+.4f\n",
  //       dr(0) - (.00509570869172363845e0),
  //       dr(1) - (.0828663025983528700e0),
  //       dr(2) - (-.0636634925404189617e0));
  
  /* Example  4 */
  printf("\nTest Case 4:\n");
  rsta.x() = 1112152.8166e0;
  rsta.y() = -4842857.5435e0;
  rsta.z() = 3985496.1783e0;
  rsun.x() = 8382471154.1312895e0;
  rsun.y() = 10512408445.356153e0;
  rsun.z() = -5360583240.3763866e0;
  rmon.x() =  380934092.93550891e0;
  rmon.y() = 2871428.1904491195e0;
  rmon.z() =  79015680.553570181e0;

  utc = TwoPartDateUTC(dso::datetime<nanoseconds>(
      year(2017), month(1), day_of_month(15), nanoseconds(0)));
  tt=(utc.utc2tt());
  ut=(tt.tt2ut1(0.5724695));

  /* Fundamental arguments (IERS 2003) */
  fa[0] = fal03(tt);
  fa[1] = falp03(tt);
  fa[2] = faf03(tt);
  fa[3] = fad03(tt);
  fa[4] = faom03(tt);
  fa[5] = fame03(tt);
  fa[6] = fave03(tt);
  fa[7] = fae03(tt);
  fa[8] = fama03(tt);
  fa[9] = faju03(tt);
  fa[10] = fasa03(tt);
  fa[11] = faur03(tt);
  fa[12] = fane03(tt);
  fa[13] = fapa03(tt);

  dr = set.displacement(tt, ut, rsta.mv, rmon.mv, rsun.mv, fa);
  //printf("Displacement : %+.6f %+.6f %+.6f\n", dr(0), dr(1), dr(2));
  //printf("Diffs        : %+.4f %+.4f %+.4f\n",
  //       dr(0) - (.0050957086917236384e0),
  //       dr(1) - (.082866302598352870e0),
  //       dr(2) - (-.063663492540418962e0 ));

  return 0;
}
