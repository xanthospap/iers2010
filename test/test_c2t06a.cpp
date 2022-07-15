#include "iau.hpp"
#include "datetime/dtfund.hpp"
#include <random>
#include <cmath>
#include <limits>

//
// test internal accuracy of c2t06a
// transform a vector using the rotation matrix produced by this function, and
// then turn the vector back and compare componants
// this program does not assert anything; it will only report results!
//

constexpr const int num_tests = 10000;

// for components of a cartesian vector
constexpr const double min_cartesian = -500e3;
constexpr const double max_cartesian = 500e3;
std::uniform_real_distribution<double> uni_car(min_cartesian, max_cartesian);

// for TT - UT1 difference
constexpr const double min_sec = 0e0;
constexpr const double max_sec = 59e0;
std::uniform_real_distribution<double> uni_sec(min_sec, max_sec);

// for xp, yp
constexpr const double min_pole = -M_PI / 5e0;
constexpr const double max_pole = M_PI / 5e0;
std::uniform_real_distribution<double> uni_pol(min_pole, max_pole);

// for MJD
constexpr const double min_mjd = 44239e0; // 01/01/1980
constexpr const double max_mjd = 59945e0; // 01/01/2023
std::uniform_real_distribution<double> uni_mjd(min_mjd, max_mjd);

int main() {
  // Seed with a real random value, if available
  std::random_device r;
  constexpr const double mind = std::numeric_limits<double>::min();
  dso::Vector3 diffs({mind, mind, mind}), avediffs({0e0, 0e0, 0e0});

  for (int t = 0; t < num_tests; t++) {
    std::default_random_engine re(r());

    // random variables ...
    const double mjd_tt = uni_mjd(re);
    const double mjd_ut1 = mjd_tt + (uni_sec(re) / 86400e0);
    const double xp = uni_pol(re);
    const double yp = uni_pol(re);

    // random vector
    dso::Vector3 vec;
    vec(0) = uni_car(re);
    vec(1) = uni_car(re);
    vec(2) = uni_car(re);
    // printf("Vector: %+22.9f %+22.9f %+22.9f\n", vec(0), vec(1), vec(2));

    // celestial - to - terrestrial matrix
    dso::Mat3x3 c2t = iers2010::sofa::c2t06a(dso::mjd0_jd, mjd_tt, dso::mjd0_jd,
                                             mjd_ut1, xp, yp);

    // apply ...
    auto tvec = c2t * vec;

    // inverse transformation
    dso::Mat3x3 t2c = c2t.transpose();
    auto cvec = t2c * tvec;

    // max-mean differences
    for (int i = 0; i < 3; i++) {
      double d;
      if ((d = std::abs(vec(i) - cvec(i))) > diffs(i)) {
        diffs(i) = d;
      }
      avediffs(i) += vec(i) - cvec(i);
    }
  }

  printf("From celestial vector to terrestrial and back. Diffs:\n");
  printf("\tDelta X Max:%.15e Average:%+.15e\n", diffs(0),
         avediffs(0) / num_tests);
  printf("\tDelta Y Max:%.15e Average:%+.15e\n", diffs(1),
         avediffs(1) / num_tests);
  printf("\tDelta Z Max:%.15e Average:%+.15e\n", diffs(2),
         avediffs(2) / num_tests);

  return 0;
}