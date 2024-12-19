#include "icgemio.hpp"
#include <charconv>
#include <cassert>

namespace {
const char *skip_ws(const char *line) noexcept {
  while (*line && *line == ' ')
    ++line;
  return line;
}

const char *skip_cols(const char *line, int num_cols) noexcept {
  const char *str = skip_ws(line);
  int cols_skipped = 0;
  while (*str && *(str + 1) && (cols_skipped < num_cols)) {
    if (*str != ' ' && *(str + 1) == ' ')
      ++cols_skipped;
    ++str;
  }
  return skip_ws(str);
}

/* translate a floating point number printed using FORTRAN's scientific 
 * notation, to a valid c++ floating point alphanumeric descriptor, i.e.
 * -.914716531750D-11 to -.914716531750e-11
 *  The 'transformed' string (at output) is stored in buf.
 *
 *  Buf should be large enough so that len characters from line are copied 
 *  to it.
 */
char *for2cpp_format(const char *line, int len, char *buf) noexcept {
  assert(len < 256);
  std::memcpy(buf, line, len * sizeof(char));
  int i = 1;
  while (i < len) {
    if ((buf[i - 1] == 'D' || buf[i - 1] == 'd') &&
        (buf[i] == '+' || buf[i] == '-')) {
      buf[i - 1] = 'e';
    }
    ++i;
  }
  return buf;
}
} /* anonymous namespace */

int dso::resolve_icgem_data_line_v1(const char *line, dso::Icgem::ErrorModel er,
                                    dso::Icgem::DataEntry &entry) noexcept {
  /* it sometimes happens (but seldom) that the line contains numeric values 
   * in FORTRAN scientific notation. Transform these (if any) to valid 
   * alphanumeric strings.
   */
  char buf[256];
  int error = 0;
  int len = std::strlen(line);
  const char *c = for2cpp_format(line, len+1, buf);

  // const char *c = line;
  int err_cols = 0;
  /* skip number of columns, depending on the Error Model of the file */
  switch (er) {
  case dso::Icgem::ErrorModel::Calibrated:
  case dso::Icgem::ErrorModel::Formal:
    err_cols = 2;
    break;
  case dso::Icgem::ErrorModel::CalibratedAndFormal:
    err_cols = 4;
    break;
  case dso::Icgem::ErrorModel::No:
    err_cols = 0;
  }

  /* Resolve Static spherical harmonic coefficients
   * We are only interested in the parameters: l,m,C,S which are always the
   * first four parameters (in this sequence).
   */
  if (!std::strncmp(line, "gfc ", 4)) {
    /* resolve degree */
    auto res = std::from_chars(skip_ws(c + 4), c + len, entry.degree);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), c + len, entry.order);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), c + len, entry.C);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), c + len, entry.S);
    if (res.ec != std::errc{})
      ++error;
    /* assign entry data */
    entry.key = Icgem::DataEntryType::gfc;

  } else if (!std::strncmp(line, "gfct", 4)) {
    /* Time variable spherical harmonic coefficients, according to the
     * respective formulae.
     *
     * degree, order, Clm, Slm, are always the first four columns;
     * t(yyyymmdd) is the last column, stored in t0.
     */
    /* resolve degree */
    auto res = std::from_chars(skip_ws(c + 4), c + len, entry.degree);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), c + len, entry.order);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), c + len, entry.C);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), c + len, entry.S);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* skip error-related cols */
    c = skip_cols(c, err_cols);
    /* resolve t at pre-last column */
    const char *ptr;
    if (resolve_icgem_data_date_v1(skip_ws(c), entry.t0, ptr))
      ++error;
    /* set key */
    entry.key = Icgem::DataEntryType::gfct;

  } else if ((!std::strncmp(line, "trnd", 4)) || (!std::strncmp(line, "dot ", 4))) {
    /* Time variable spherical harmonic coefficients, according to the
     * respective formulae.
     *
     * degree, order, Clm, Slm, are always the first four columns;
     */
    /* resolve degree */
    auto res = std::from_chars(skip_ws(c + 4), c + len, entry.degree);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), c + len, entry.order);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), c + len, entry.C);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), c + len, entry.S);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* set key */
    entry.key = Icgem::DataEntryType::trnd;

  } else if (!std::strncmp(line, "asin", 4)) {
    /* Time variable spherical harmonic coefficients, according to the
     * respective formulae.
     *
     * degree, order, sine_amplitude_C, sine_amplitude_S,
     * are always the first four columns;
     *
     * last column is period [year]
     */
    /* resolve degree */
    auto res = std::from_chars(skip_ws(c + 4), c + len, entry.degree);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), c + len, entry.order);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), c + len, entry.C);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), c + len, entry.S);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* skip error-related cols */
    c = skip_cols(c, err_cols);
    /* resolve period */
    res = std::from_chars(skip_ws(c), c + len, entry.period);
    if (res.ec != std::errc{})
      ++error;
    /* set key */
    entry.key = Icgem::DataEntryType::asin;

  } else if (!std::strncmp(line, "acos", 4)) {
    /* Time variable spherical harmonic coefficients, according to the
     * respective formulae.
     *
     * degree, order, sine_amplitude_C, sine_amplitude_S,
     * are always the first four columns;
     *
     * last column is period [year]
     */
    /* resolve degree */
    auto res = std::from_chars(skip_ws(c + 4), c + len, entry.degree);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), c + len, entry.order);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), c + len, entry.C);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), c + len, entry.S);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* skip error-related cols */
    c = skip_cols(c, err_cols);
    /* resolve period */
    res = std::from_chars(skip_ws(c), c + len, entry.period);
    if (res.ec != std::errc{})
      ++error;
    entry.key = Icgem::DataEntryType::acos;

  } else {
    fprintf(stderr,
            "[ERROR] Unknown keyword encountered in the ICGEM file! line [%s], "
            "(traceback: %s)\n",
            line, __func__);
    ++error;
  }

  return error;
}
