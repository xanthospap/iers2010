#include "fcn.hpp"
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <datetime/tpdate.hpp>
#include <fstream>

/** Structure of a coefficient file:
 * %           mjd     real     imag    sigma  (all in microas)
 * 1984.0  45708.0   -220.93   -147.19     10.22
 * 1984.5  45890.6   -219.08   -144.10      9.76
 * 1985.0  46073.2   -208.97   -103.48      9.29
 * 1985.5  46255.9   -213.49   -110.42      8.63
 * 1986.0  46438.5   -185.10   -111.47      8.10
 * 1986.5  46621.1   -197.09   -108.99      7.48
 * 1987.0  46803.8   -174.08    -76.14      6.75
 * 1987.5  46986.4   -187.09    -53.60      6.29
 * 1988.0  47169.0   -174.22    -57.70      5.79
 */

namespace {
constexpr const std::size_t LINE_LEN = 128;
inline const char *skip_ws(const char *line) noexcept {
  while (*line && *line == ' ')
    ++line;
  return line;
}
} /* unnamed namespace */

int dso::parse_lambert_coefficients(
    const char *fn, const dso::TwoPartDate &from, const dso::TwoPartDate &to,
    std::vector<fcn::LambertCoefficients> &lvec) noexcept {
  if (to < from) {
    fprintf(stderr,
            "[ERROR] Invalid time interval request for parsing Lambert "
            "coefficients (traceback: %s)\n",
            __func__);
    return 5;
  }
  
  /* clear and allocated space for the coefficients entry */
  lvec.clear();
  /* approximatelly two entries per year */
  lvec.reserve(std::ceil(
      (to.diff<dso::DateTimeDifferenceType::FractionalYears>(from) + 1) * 2e0));

  /* open the coefficients file */
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(
        stderr,
        "[ERROR] Failed to open Lambert coefficients file %s (traceback: %s)\n",
        fn, __func__);
    return 1;
  }

  char line[LINE_LEN];

  /* read and discard the first line; should be the header starting with a
   * '%' character
   */
  if (!fin.getline(line, LINE_LEN)) {
    fprintf(stderr,
            "[ERROR] Failed reading Lambert coefficients header from file %s "
            "(traceback: %s)\n",
            fn, __func__);
    return 1;
  }
  if (*skip_ws(line) != '%') {
    fprintf(stderr,
            "[ERROR] Failed reading Lambert coefficients header from file %s; "
            "expected character \'%%\' (traceback: %s)\n",
            fn, __func__);
    return 2;
  }

  int error = 0;
  const int max_allowed_lines = 1000;
  int line_nr = 0;
  double d[5];
  /* loop through file lines */
  while ((fin.getline(line, LINE_LEN)) && (!error) &&
         (++line_nr < max_allowed_lines)) {
    const char *str = line;
    int sz = std::strlen(str);
    /* resolve the numeric values */
    for (int i = 0; i < 5; i++) {
      const auto res = std::from_chars(skip_ws(str), line + sz, d[i]);
      if (res.ec != std::errc())
        ++error;
      str = res.ptr;
    }
    /* are we interested ? */
    if (d[1] >= from.as_mjd() && d[1] < to.as_mjd()) {
      /* resolve fractional MJD to a TwoPartDate */
      double ipart;
      const double fday = std::modf(d[1], &ipart);
      const dso::TwoPartDate t((int)ipart, fday * dso::SEC_PER_DAY);
      /* push back the coefficients entry */
      lvec.emplace_back(t, d[2], d[3], d[4]);
    } else {
      ;
    }
  }

  /* check for errors */
  if (error || (line_nr >= max_allowed_lines)) {
    fprintf(stderr,
            "[ERROR] Failed reading Lambert coefficients from file %s "
            "(traceback: %s)\n",
            fn, __func__);
    return 3;
  }
  if (!fin && !fin.eof()) {
    fprintf(stderr,
            "[ERROR] Failed reading Lambert coefficients from file %s "
            "(traceback: %s)\n",
            fn, __func__);
    return 4;
  }

  /* all done */
  return 0;
}
