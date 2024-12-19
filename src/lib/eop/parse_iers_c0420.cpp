#include "eop.hpp"
#include <charconv>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <cassert>

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

/*
 * # EARTH ORIENTATION PARAMETER (EOP) PRODUCT CENTER CENTER (PARIS OBSERVATORY)
 * - INTERNATIONAL EARTH ROTATION AND REFERENCE SYSTEMS SERVICE # EOP (IERS) 20
 * C04 TIME SERIES  consistent with ITRF 2020 - sampled at 0h UTC # Contact:
 * christian.bizouard@obspm.fr # Reference Precession-Nutation Model: IAU 2000 #
 * format(4(i4),f10.2,2(f12.6),f12.7,2(f12.6),2(f12.6),f12.7,2(f12.6),f12.7,2(f12.6),2(f12.6),f12.7)
 * # YR  MM  DD  HH       MJD        x(")        y(")  UT1-UTC(s)       dX(")
 * dY(") xrt(")      yrt(")      LOD(s)        x Er        y Er  UT1-UTC Er dX
 * Er dY Er       xrt Er      yrt Er      LOD Er 1962   1   1   0  37665.00
 * -0.012700 0.213000   0.0326338    0.000000    0.000000    0.000000 0.000000
 * 0.0017230 0.030000    0.030000   0.0020000    0.004774    0.002000 0.000000
 * 0.000000 0.0014000 1962   1   2   0  37666.00   -0.015900    0.214100
 * 0.0320547 0.000000    0.000000    0.000000    0.000000   0.0016690 0.030000
 * 0.030000 0.0020000    0.004774    0.002000    0.000000    0.000000 0.0014000
 */
int dso::details::parse_iers_C0420(const char *c04fn,
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
    int days = std::ceil(
        end.diff<dso::DateTimeDifferenceType::FractionalDays>(start).days());
    eops.clear();
    eops.reserve(days);
  }

  char line[MAX_LINE_CHARS];
  dso::EopRecord rec;
  int error = 0;
  /* keep on reading/parsing as long as we don;t encounter a date > end */
  while (fin.getline(line, MAX_LINE_CHARS) && !error) {
    const char *cstart;
    const char *cend = line + std::strlen(line);
    if (line[0] && (*line != ' ' && *line != '#')) {
      /* skip the data (YYYY MM DD) which is 12 chars length and get the mjd
       * (note that this is UTC)
       */
      cstart = next_num(line + 12);
      std::from_chars_result fcr;

      /* parse date -- made up or recorded hours (in day) and MJD */
      dso::TwoPartDateUTC ct;
      { /* we have to read in (integer) hours */
        int hours;
        fcr = std::from_chars(cstart, cend, hours);
        if (fcr.ec != std::errc() || fcr.ptr == cstart) {
          fprintf(
              stderr,
              "[ERROR] Failed reading hours from eop file %s (traceback: %s)\n",
              c04fn, __func__);
          return 2;
        }
        cstart = next_num(fcr.ptr);
        /* read real MJD */
        double mjd;
        fcr = std::from_chars(cstart, cend, mjd);
        if (fcr.ec != std::errc() || fcr.ptr == cstart) {
          fprintf(
              stderr,
              "[ERROR] Failed reading MJD from eop file %s (traceback: %s)\n",
              c04fn, __func__);
          return 2;
        }
        cstart = next_num(fcr.ptr);
        /* split the real MJD to integral and fractional part, and verify */
        double imjd;
        double fdays = std::modf(mjd, &imjd);
        assert(fdays == hours/24e0);
        ct = dso::TwoPartDateUTC(
            (int)imjd,
            dso::FractionalSeconds{fdays * dso::SEC_PER_DAY});
      }

      /* UTC date to TT and ΔAT at this UTC */
      const dso::TwoPartDate ctt = ct.utc2tt();
      const auto dat = dso::dat(dso::modified_julian_day(ct.imjd()));
      /* if the date is ok, collect data */
      if (ctt >= start && ctt < end) {
        rec.t() = ctt;
        /* collecting the following (values/units)
         *  x  y   UT1-UTC   dX   dY  xrt     yrt     LOD
         *  "  "   sec       "    "   "/day   "/day   sec
         */
        double data[8];
        for (int i = 0; i < 8; i++) {
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
        rec.dX() = data[3];
        rec.dY() = data[4];
        rec.xp_rate() = data[5];
        rec.yp_rate() = data[6];
        rec.lod() = data[7];
        rec.dat() = dat;

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
