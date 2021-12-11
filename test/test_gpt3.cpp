#include "tropo.hpp"
#include <cmath>
#include <cstdio>
#include <random>

double deg2rad(double rad) noexcept { return rad * M_PI / 180e0; }
double rad2deg(double deg) noexcept { return deg * 180e0 / M_PI; }

struct vmf3_mf {
  double mfh, mfw;
};

std::random_device rd; // only used once to initialise (seed) engine
std::mt19937
    rng(rd()); // random-number engine used (Mersenne-Twister in this case)
std::uniform_int_distribution<int> uni(0e3, 3e6); // guaranteed unbiased
std::uniform_real_distribution<double> runi(0e3,
                                            M_PI / 2e0); // guaranteed unbiased

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "ERROR! Usage: %s [GPT3-NxN GRID FILE]\n", argv[0]);
    return 1;
  }

  long num_sta = 0;
  double lat0 = -90e0;
  double lat1 = 90e0;
  double lon0 = -180e0;
  double lon1 = 360e0;
  double lat_step = 3.33;
  double lon_step = 4.33;

  for (double lat = lat0; lat <= lat1; lat += lat_step) {
    for (double lon = lon0; lon <= lon1; lon += lon_step) {
      ++num_sta;
    }
  }

  double *latv = new double[num_sta];
  double *lonv = new double[num_sta];
  double *hell = new double[num_sta];
  double *znth = new double[num_sta];
  vmf3_mf *vmf3 = new vmf3_mf[num_sta];
  dso::vmf3_hw *vhw = new dso::vmf3_hw[num_sta];
  dso::gpt3_result *res = new dso::gpt3_result[num_sta];

  long ix = 0;
  for (double lat = lat0; lat <= lat1; lat += lat_step) {
    for (double lon = lon0; lon <= lon1; lon += lon_step) {
      latv[ix] = deg2rad(lat);
      lonv[ix] = deg2rad(lon);
      hell[ix] = uni(rng) / 1e3;
      znth[ix++] = runi(rng);
    }
  }

  dso::datetime<dso::nanoseconds> t{
      dso::year{2021}, dso::month{11}, dso::day_of_month{12},
      dso::nanoseconds{(6 * 60 * 60L + 123456L) *
                       dso::nanoseconds::sec_factor<long>()}};

  int grid_step = 0;
  const char *grid = argv[1];
  int error =
      dso::gpt3_fast(t, latv, lonv, hell, num_sta, 1, grid, res, grid_step);

  if (error)
    return error;

  double mfh, mfw;
  for (int i = 0; i < num_sta; i++) {
    printf("%.12f %+.12f %+.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f "
           "%.12f %.12f %.12f %.12f "
           "%.12f %.12f\n",
           t.as_mjd(), rad2deg(latv[i]), rad2deg(lonv[i]), res[i].p, res[i].T,
           res[i].dT, res[i].Tm, res[i].e, res[i].ah, res[i].aw, res[i].la,
           res[i].undu, res[i].Gn_h, res[i].Ge_h, res[i].Gn_w, res[i].Ge_w);
    dso::vmf3(res[i].ah, res[i].aw, t, latv[i], lonv[i], znth[i], mfh, mfw);
    vmf3[i] = vmf3_mf{mfh, mfw};
  }

  // print vmf3 results in a seperate file
  FILE *fptr = fopen("vmf3_res.cpp", "w");
  for (int i = 0; i < num_sta; i++) {
    fprintf(fptr, "%.12f %.12f\n", vmf3[i].mfh, vmf3[i].mfw);
  }
  fclose(fptr);

  // check results between one-station-ata-time and vectorized call
  dso::vmf3(res, t, latv, lonv, znth, vhw, num_sta);
  for (int i = 0; i < num_sta; i++) {
    double diff = std::abs(vmf3[i].mfh - vhw[i].mfh);
    if (diff > 1e-15) {
      fprintf(stderr, ">> discrepancy >1e15 for hydrostatic component!\n");
    }
    diff = std::abs(vmf3[i].mfw - vhw[i].mfw);
    if (diff > 1e-15) {
      fprintf(stderr, ">> discrepancy >1e15 for wet component!\n");
    }
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
  fprintf(fp, "grid=gpt3_%d_fast_readGrid();\n", grid_step);
  fprintf(fp,
          "[p,T,dT,Tm,e,ah,aw,la,undu,Gn_h,Ge_h,Gn_w,Ge_w] = gpt3_%d_fast(mjd, "
          "lat, lon, h_ell, %d, grid);\n",
          grid_step, 1);
  fprintf(fp, "fid = fopen(\"octave_gpt3_results.txt\", \"w\");\n");
  fprintf(fp, "for i=1:%d\n", (int)num_sta);
  fprintf(fp,
          "fprintf(fid, \"%%.12f %%+.12f %%+.12f %%.12f %%.12f %%.12f "
          "%%.12f %%.12f %%.12f %%.12f %%.12f %%.12f %%.12f %%.12f "
          "%%.12f %%.12f\\n\", mjd, rad2deg(lat(i)), rad2deg(lon(i)), p(i), "
          "T(i), dT(i), Tm(i), e(i), ah(i), aw(i), la(i), undu(i), Gn_h(i), "
          "Ge_h(i), Gn_w(i), Ge_w(i));\n");
  fprintf(fp, "endfor\n");
  fprintf(fp, "fclose(fid)\n");
  fprintf(fp, "vmfh=zeros(%d,1); vmfw=zeros(%d,1);\n", (int)num_sta,
          (int)num_sta);
  for (int i = 0; i < num_sta; i++) {
    fprintf(fp,
            "[vmfh(%d), vmfw(%d)] = vmf3(%.20f, %.20f, %.20f, %20f, %.20f, "
            "%.20f);\n",
            i + 1, i + 1, res[i].ah, res[i].aw, t.as_mjd(), latv[i], lonv[i],
            znth[i]);
  }
  fprintf(fp, "fid = fopen(\"octave_vmf3_results.txt\", \"w\");\n");
  fprintf(fp, "for i=1:%d\n", (int)num_sta);
  fprintf(fp, "\tfprintf(fid, \"%%.15f %%.15f\\n\", vmfh(i), vmfw(i));");
  fprintf(fp, "end\n");
  fprintf(fp, "fclose(fid)\n");
  fclose(fp);

  printf("Note that the input grid file has a step of %d degrees\n", grid_step);
  return error;
}
