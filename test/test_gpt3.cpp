#include "gpt3.hpp"
#include <random>

std::random_device rd; // only used once to initialise (seed) engine
std::mt19937
    rng(rd()); // random-number engine used (Mersenne-Twister in this case)
std::uniform_int_distribution<int> uni(0e3, 3e6); // guaranteed unbiased

int main() {
  long num_sta = 0;

  for (double lat = -M_PI / 2e0; lat <= M_PI / 2e0; lat += 0.333333) {
    for (double lon = -M_PI; lon <= 2e0 * M_PI; lon += 0.533533) {
      ++num_sta;
    }
  }

  double *latv = new double[num_sta];
  double *lonv = new double[num_sta];
  double *hell = new double[num_sta];
  dso::gpt3_result *res = new dso::gpt3_result[num_sta];

  long ix = 0;
  for (double lat = -M_PI / 2e0; lat <= M_PI / 2e0; lat += 0.333333) {
    for (double lon = -M_PI; lon <= 2e0 * M_PI; lon += 0.533533) {
      latv[ix] = lat;
      lonv[ix] = lon;
      hell[ix] = uni(rng) / 1e3;
    }
  }

  dso::datetime<dso::nanoseconds> t{
      dso::year{2021}, dso::month{11}, dso::day_of_month{12},
      dso::nanoseconds{(6 * 60 * 60L + 123456L) *
                       dso::nanoseconds::sec_factor<long>()}};

  const char *grid =
      "/home/xanthos/Software/iers2010/test/gpt3/gpt3_5.grd";
  int error =
          gpt3_5_fast(t, latv, lonv, hell, num_sta, 0, grid, res);

  for (int i = 0; i < num_sta; i++) {
    printf("%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f "
           "%.12f %.12f\n",
           t.as_mjd(), res[i].p, res[i].T, res[i].dT, res[i].Tm, res[i].e,
           res[i].aw, res[i].ah, res[i].undu, res[i].Gn_h,
           res[i].Ge_h, res[i].Gn_w, res[i].Gn_w);
  }

  return error;
}