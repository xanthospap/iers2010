#include "gpt3.hpp"
#include <cmath>
#include <random>
#include <cstdio>

double deg2rad(double rad) noexcept { return rad * M_PI / 180e0; }
double rad2deg(double deg) noexcept { return deg * 180e0 / M_PI; }

std::random_device rd; // only used once to initialise (seed) engine
std::mt19937
    rng(rd()); // random-number engine used (Mersenne-Twister in this case)
std::uniform_int_distribution<int> uni(0e3, 3e6); // guaranteed unbiased

int main() {
  long num_sta = 0;
  double lat0 = -90e0;
  double lat1 = 90e0;
  double lon0 = -180e0;
  double lon1 = 360e0;
  double lat_step = 30.33;
  double lon_step = 46.33;

  for (double lat = lat0; lat <= lat1; lat += lat_step) {
    for (double lon = lon0; lon <= lon1; lon += lon_step) {
      ++num_sta;
    }
  }

  double *latv = new double[num_sta];
  double *lonv = new double[num_sta];
  double *hell = new double[num_sta];
  dso::gpt3_result *res = new dso::gpt3_result[num_sta];

  long ix = 0;
  for (double lat = lat0; lat <= lat1; lat += lat_step) {
    for (double lon = lon0; lon <= lon1; lon += lon_step) {
      latv[ix] = deg2rad(lat);
      lonv[ix] = deg2rad(lon);
      hell[ix++] = uni(rng) / 1e3;
    }
  }

  dso::datetime<dso::nanoseconds> t{
      dso::year{2021}, dso::month{11}, dso::day_of_month{12},
      dso::nanoseconds{(6 * 60 * 60L + 123456L) *
                       dso::nanoseconds::sec_factor<long>()}};

  const char *grid =
      //      "/home/xanthos/Software/iers2010/test/gpt3/gpt3_5.grd";
      "/home/xanthos/Builds/iers2010/test/gpt3_5.grd";
  int error = gpt3_5_fast(t, latv, lonv, hell, num_sta, 1, grid, res);

  if (error)
    return error;

  for (int i = 0; i < num_sta; i++) {
    printf("%.12f %+.12f %+.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f "
           "%.12f %.12f %.12f %.12f "
           "%.12f %.12f\n",
           t.as_mjd(), rad2deg(latv[i]), rad2deg(lonv[i]), res[i].p, res[i].T,
           res[i].dT, res[i].Tm, res[i].e, res[i].ah, res[i].aw, res[i].la, res[i].undu,
           res[i].Gn_h, res[i].Ge_h, res[i].Gn_w, res[i].Ge_w);
  }

  // write an octave script to check results
  FILE *fp = fopen("test_gpt3.m", "w");
  fprintf(fp, "%% Octave script %%\n");
  fprintf(fp, "lat=[");
  for (int i = 0; i < num_sta; i++)
    fprintf(fp, "%.15f;", latv[i]);
  fprintf(fp, "];\n");
  fprintf(fp, "lon=[");
  for (int i = 0; i < num_sta; i++)
    fprintf(fp, "%.15f;", lonv[i]);
  fprintf(fp, "];\n");
  fprintf(fp, "h_ell=[");
  for (int i = 0; i < num_sta; i++)
    fprintf(fp, "%.15f;", hell[i]);
  fprintf(fp, "];\n");
  fprintf(fp, "mjd=%.15f\n", t.as_mjd());
  fprintf(fp, "grid=gpt3_5_fast_readGrid();\n");
  fprintf(fp, "[p,T,dT,Tm,e,ah,aw,la,undu,Gn_h,Ge_h,Gn_w,Ge_w] = gpt3_5_fast(mjd, "
         "lat, lon, h_ell, %d, grid);\n",
         1);
  fprintf(fp, "fid = fopen(\"octave_gpt3_results.txt\", \"w\");\n");
      fprintf(fp, "for i=1:%d\n", (int)num_sta);
      fprintf(fp, "fprintf(fid, \"%%.12f %%+.12f %%+.12f %%.12f %%.12f %%.12f "
                  "%%.12f %%.12f %%.12f %%.12f %%.12f %%.12f %%.12f %%.12f "
                  "%%.12f %%.12f\\n\", mjd, rad2deg(lat(i)), rad2deg(lon(i)), p(i), "
                  "T(i), dT(i), Tm(i), e(i), ah(i), aw(i), la(i), undu(i), Gn_h(i), "
                  "Ge_h(i), Gn_w(i), Ge_w(i));\n");
      fprintf(fp, "endfor\n");
      fprintf(fp, "fclose(fid)\n");
      fclose(fp);

      return error;
}
