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

  double dr, dw, ds;
  dso::Hardisp h;
  h(mjd);
  h.hardisp(blqInfoVec[0], dr, dw, ds);
  
  printf("dUp: %+.6f dSouth: %+.6f dWest: %+.6f\n",dr, ds, dw);

  // compute displacement for the next six hours
  for (int i = 0; i < 24; i++) {
    t.add_seconds(
        dso::cast_to<dso::seconds, dso::nanoseconds>(dso::seconds(3600L)));
    dso::TwoPartDate d(t);
    h(d);
    h.hardisp(blqInfoVec[0], dr, dw, ds);
    printf("dUp: %+.6f dSouth: %+.6f dWest: %+.6f\n",dr, ds, dw);
  }

  return 0;
}
