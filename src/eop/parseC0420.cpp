#include "eop.hpp"
#include <charconv>
#include <cstdio>
#include <cstring>
#include <fstream>

namespace {
constexpr const std::size_t MAX_LINE_CHARS = 256;

/// @brief Skip whitespace characters. Return the first substring the does not
///        start with a ' ' character.
const char* next_num(const char* line) noexcept
{
  const char* ptr = line;
  while (*ptr && *ptr == ' ')
    ++ptr;
  return ptr;
}
} // unnamed namespace

/*
# EARTH ORIENTATION PARAMETER (EOP) PRODUCT CENTER CENTER (PARIS OBSERVATORY) -
INTERNATIONAL EARTH ROTATION AND REFERENCE SYSTEMS SERVICE # EOP (IERS) 20 C04
TIME SERIES  consistent with ITRF 2020 - sampled at 0h UTC # Contact:
christian.bizouard@obspm.fr # Reference Precession-Nutation Model: IAU 2000 #
format(4(i4),f10.2,2(f12.6),f12.7,2(f12.6),2(f12.6),f12.7,2(f12.6),f12.7,2(f12.6),2(f12.6),f12.7)
# YR  MM  DD  HH       MJD        x(")        y(")  UT1-UTC(s)       dX(") dY(")
xrt(")      yrt(")      LOD(s)        x Er        y Er  UT1-UTC Er      dX Er dY
Er       xrt Er      yrt Er      LOD Er 1962   1   1   0  37665.00   -0.012700
0.213000   0.0326338    0.000000    0.000000    0.000000    0.000000   0.0017230
0.030000    0.030000   0.0020000    0.004774    0.002000    0.000000    0.000000
0.0014000 1962   1   2   0  37666.00   -0.015900    0.214100   0.0320547
0.000000    0.000000    0.000000    0.000000   0.0016690    0.030000    0.030000
0.0020000    0.004774    0.002000    0.000000    0.000000   0.0014000
*/
int dso::parse_iers_C0420(const char* c04fn, dso::modified_julian_day start_mjd,
    dso::modified_julian_day end_mjd,
    dso::EopLookUpTable& eoptable) noexcept
{
  // open file
  std::ifstream fin(c04fn);
  if (!fin.is_open()) {
    fprintf(stderr,
        "[ERROR] Failed opening EOP (C04) file %s (traceback: "
        "%s)\n",
        c04fn, __func__);
    return 1;
  }

  // if needed, resize the EOP table, interval is [start.end)
  int days = end_mjd.as_underlying_type() - start_mjd.as_underlying_type();
  eoptable.clear();
  eoptable.reserve(days);

  // start and end MJD as dates
  const dso::TwoPartDate start_t(
      static_cast<double>(start_mjd.as_underlying_type()), 0e0);
  const dso::TwoPartDate end_t(
      static_cast<double>(end_mjd.as_underlying_type()), 0e0);

  char line[MAX_LINE_CHARS];

  dso::EopRecord rec;
  while (fin.getline(line, MAX_LINE_CHARS)) {
    const char* start;
    const char* end = line + std::strlen(line);
    if (*line && (*line != ' ' && *line != '#')) {
      // skip the data (YYYY MM DD) which is 12 chars length and get the mjd
      // (note that this is UTC)
      start = next_num(line + 12);
      std::from_chars_result fcr;

      // parse date -- made up or recorded hours (in day) and MJD
      dso::TwoPartDate ct;
      { // we have to read in (integer) hours
        int hours;
        fcr = std::from_chars(start, end, hours);
        if (fcr.ec != std::errc() || fcr.ptr == start) {
          fprintf(
              stderr,
              "[ERROR] Failed reading hours from eop file %s (traceback: %s)\n",
              c04fn, __func__);
          return 2;
        }
        start = next_num(fcr.ptr);
        const double fractional_days = hours / 24e0;
        double mjd;
        fcr = std::from_chars(start, end, mjd);
        if (fcr.ec != std::errc() || fcr.ptr == start) {
          fprintf(
              stderr,
              "[ERROR] Failed reading MJD from eop file %s (traceback: %s)\n",
              c04fn, __func__);
          return 2;
        }
        start = next_num(fcr.ptr);
        ct = dso::TwoPartDate(mjd, fractional_days).normalized();
      }

      int error = 0;

      // if the date is ok, collect data
      if (ct >= start_t && ct < end_t) {
        error = 0;
        rec.mjd = ct;
        //  x  y   UT1-UTC   dX   dY  xrt     yrt     LOD
        //  "  "   sec       "    "   "/day   "/day   sec
        double data[8];
        for (int i = 0; i < 8; i++) {
          start = next_num(fcr.ptr);
          fcr = std::from_chars(start, end, data[i]);
          if (fcr.ec != std::errc() || fcr.ptr == start) {
            ++error;
          }
        }

        if (error) {
          fprintf(stderr,
              "[ERROR] Failed resolving EOP line [%s]. EOP file is %s "
              "(traceback: %s/err=B)\n",
              line, c04fn, __func__);
          return 2;
        }

        // assign (watch the order)
        rec.xp = data[0];
        rec.yp = data[1];
        rec.dut = data[2];
        rec.dx = data[3];
        rec.dy = data[4];
        rec.xrt = data[5];
        rec.yrt = data[6];
        rec.lod = data[7];

        eoptable.push_back(rec);
      } else if (ct >= end_t) {
        break;
      }

    } // if (line[0] != ' ') ...
  }   // end reading lines

  // we were supposed to collect:
  if (days != eoptable.size()) {
    fprintf(stderr,
        "[WRNNG] Failed to collect EOP data for requested date span\n");
    if (!eoptable.size()) {
      fprintf(stderr, "[WRNNG] No epoch collected from EOP file!\n");
    } else {
      fprintf(stderr,
          "[WRNNG] Requested [%ld to %ld) collected [%.3f to %.3f]\n",
          start_mjd.as_underlying_type(), end_mjd.as_underlying_type(),
          eoptable.epoch_vector()[0].mjd(),
          eoptable.epoch_vector().back().mjd());
    }
  }
  return !(days == eoptable.size());
}
