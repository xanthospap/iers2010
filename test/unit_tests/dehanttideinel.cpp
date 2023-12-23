#include "solid_earth_tide.hpp"
#include "fundarg.hpp"
#include <cstdio>

using namespace dso;

int main() {

  dso::CartesianCrd rsta, rmon, rsun;
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
  const double fa[14] = {
      fal03(tt),  falp03(tt),
      faf03(tt),  fad03(tt),
      faom03(tt), fame03(tt),
      fave03(tt), fae03(tt),
      fama03(tt), faju03(tt),
      fasa03(tt), faur03(tt),
      fane03(tt), fapa03(tt),
  };

  SolidEarthTide set(3.986004418e14, 6378136.6e0, 4.9048695e12,
                     1.32712442099e20);

  const auto dr = set.displacement(tt, ut, rsta.mv, rmon.mv, rsun.mv, fa);

  printf("Displacement coefficients: %+.6f %+.6f %+.6f\n", dr(0), dr(1), dr(2));

  return 0;
}
