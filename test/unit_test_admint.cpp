#include "hardisp.hpp"
#include <cstdio>
#include "datetime/dtcalendar.hpp"
#include <vector>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [BLQ FILE]\n", argv[0]);
    return 1;
  }

  std::vector<dso::BlqSiteInfo> blqInfoVec;
  std::vector<const char *> sites;
  const char *site = "ONSA";
  sites.emplace_back(site);
  if (parse_blq(argv[1], blqInfoVec, &sites)) {
    fprintf(stderr, "Failed reading BLQ file\n");
    return 1;
  }

  // 2009 6 25 1 10 45
  dso::datetime<dso::nanoseconds> t(
      dso::year(2009), dso::month(6), dso::day_of_month(25), dso::hours(1),
      dso::minutes(10), dso::nanoseconds(45L * 1'000'000'000L));
  dso::TwoPartDate mjd(t);

  // dso::AdmintReturnStruct res;
  dso::Hardisp h;
  h(mjd);
  h.admint(blqInfoVec[0]);

  const auto res = h.admitances();
  for (int i = 0; i < (int)res.vadmin.size(); i++) {
    printf("%3d F: %.6f AMP: %9.6f%9.6f%9.6f PH: %12.6f%12.6f%12.6f\n", i,
           res.vadmin[i].freq,
           res.vadmin[i].ampl_r, res.vadmin[i].ampl_w, res.vadmin[i].ampl_s,
           dso::rad2deg(res.vadmin[i].ph_r), dso::rad2deg(res.vadmin[i].ph_w),
           dso::rad2deg(res.vadmin[i].ph_s));
  }

  return 0;
}
