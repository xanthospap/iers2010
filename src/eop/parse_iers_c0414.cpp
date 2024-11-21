#include "eop.hpp"
#include <charconv>
#include <cstdio>
#include <cstring>
#include <fstream>

namespace {
constexpr const std::size_t MAX_LINE_CHARS = 256;

/* @brief Skip whitespace characters. Return the first substring the does not
 *        start with a ' ' character.
 */
const char *next_num(const char *line) noexcept {
  const char *ptr = line;
  while (ptr && *ptr == ' ')
    ++ptr;
  return ptr;
}
} /* unnamed namespace */

/* Example : https://datacenter.iers.org/data/224/eopc04_14_IAU2000.62-now.txt
 * Format of data:
 * FORMAT(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2(F12.6))
  Date      MJD      x          y        UT1-UTC       LOD         dX dY x Err y Err   UT1-UTC Err  LOD Err     dX Err       dY Err "          "           s s
"         "           "          "          s         s            "           "
     (0h UTC)

1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230   0.000000
0.000000   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000
1962   1   2  37666  -0.015900   0.214100   0.0320547   0.0016690   0.000000
0.000000   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000
*/
int dso::details::parse_iers_C0414(const char *c04fn,
                                   const dso::MjdEpoch &start,
                                   const dso::MjdEpoch &end,
                                   dso::EopSeries &eops) noexcept {
  /* check dates */
  if (end <= start) {
    fprintf(stderr,
            "[ERROR] Invalid dates for parsing EOP (C04) file %s (traceback: "
            "%s)\n",
            c04fn, __func__);
    return 1;
  }

  /* open file */
  std::ifstream fin(c04fn);
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed opening EOP (C04) file %s (traceback: "
            "%s)\n",
            c04fn, __func__);
    return 1;
  }

  {
    /* if needed, resize the EOP table, interval is [start, end) */
    int days =
        std::ceil(end.diff<dso::DateTimeDifferenceType::FractionalDays>(start).days());
    eops.clear();
    eops.reserve(days);
  }

  char line[MAX_LINE_CHARS];
  dso::EopRecord rec;
  int error=0;
  /* keep on reading lines untill we encounter a date later than end */
  while (fin.getline(line, MAX_LINE_CHARS) && !error) {
    const char *cstart;
    const char *cend = line + std::strlen(line);
    long imjd=0;
    if (line[0] && (*line != ' ' && *line != '#')) {
      /* skip the data (YYYY MM DD) which is 12 chars length and get the mjd
       * (note that this is UTC)
       */
      cstart = next_num(line + 12);
      auto fcr = std::from_chars(cstart, cend, imjd);
      if (fcr.ec != std::errc() || fcr.ptr == cstart)
        error = 1;

      /* integral part of current MJD (UTC) */
      const dso::TwoPartDateUTC cutc(imjd, dso::FractionalSeconds{0e0});
      const dso::TwoPartDate ctt = cutc.utc2tt();
      /* ΔAT at time cutc */
      const auto dat = dso::dat(dso::modified_julian_day(cutc.imjd()));

      /* if the date is within range, collect data */
      if (ctt >= start && ctt < end) {
        rec.t() = ctt;
        /* Collect the following (in given order/units :
         * x, y, UT1-UTC, LOD, dX, dY, [..errors...]
         * "  "  sec      sec  "   "
         */
        double data[6];
        for (int i = 0; i < 6; i++) {
          cstart = next_num(fcr.ptr);
          fcr = std::from_chars(cstart, cend, data[i]);
          if (fcr.ec != std::errc() || fcr.ptr == cstart) {
            ++error;
          }
        }

        /* assign (watch the order) */
        rec.xp() = data[0];
        rec.yp() = data[1];
        rec.dut() = data[2];
        rec.lod() = data[3];
        rec.dX() = data[4];
        rec.dY() = data[5];
        rec.dat() = dat;
        /* no x/y rate availabel, set to 0 */
        rec.xp_rate() = 0e0;
        rec.yp_rate() = 0e0;

        /* add entries to the EopLookUpTable */
        if (!error)
          eops.push_back(rec);
      } else if (ctt >= end) {
        break;
      }

    } /* if (line[0] != ' ') ... */
  }   /* end reading lines */

  if (error || (!fin && !fin.eof())) {
    fprintf(stderr,
            "[ERROR] Failed resolving EOP line [%s]. EOP file is %s "
            "(traceback: %s)\n",
            line, c04fn, __func__);
    return 1;
  }

  /* all done */
  return 0;
}
