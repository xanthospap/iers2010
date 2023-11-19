#include "icgemio.hpp"
#include <charconv>

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
} /* anonymous namespace */

int dso::resolve_icgem_data_line_v2(const char *line, dso::Icgem::ErrorModel er,
                                    dso::Icgem::DataEntry &entry) noexcept {
  int error = 0;
  int len = std::strlen(line);
  const char *c = line;
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
    auto res = std::from_chars(skip_ws(c + 4), line + len, entry.degree);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), line + len, entry.order);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), line + len, entry.C);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), line + len, entry.S);
    if (res.ec != std::errc{})
      ++error;
    /* assign entry data */
    entry.key = Icgem::DataEntryType::gfc;

  } else if (!std::strncmp(line, "gfct", 4)) {
    /* Time variable spherical harmonic coefficients, according to the
     * respective formulae.
     *
     * degree, order, Clm, Slm, are always the first four columns;
     * t0(yyyymmdd.hhmm) and t1(yyyymmdd.hhmm) are always the last two
     * columns.
     */
    /* resolve degree */
    auto res = std::from_chars(skip_ws(c + 4), line + len, entry.degree);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), line + len, entry.order);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), line + len, entry.C);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), line + len, entry.S);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* skip error-related cols */
    c = skip_cols(c, err_cols);
    /* resolve t0 at pre-last column */
    const char *ptr;
    if (resolve_icgem_data_date_v2(skip_ws(c), entry.t0, ptr))
      ++error;
    /* resolve t1 at last column */
    if (!error)
      error += resolve_icgem_data_date_v2(skip_ws(ptr), entry.t1, ptr);
    /* set key */
    entry.key = Icgem::DataEntryType::gfct;

  } else if (!std::strncmp(line, "trnd", 4)) {
    /* Time variable spherical harmonic coefficients, according to the
     * respective formulae.
     *
     * degree, order, Clm, Slm, are always the first four columns;
     */
    /* resolve degree */
    auto res = std::from_chars(skip_ws(c + 4), line + len, entry.degree);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), line + len, entry.order);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), line + len, entry.C);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), line + len, entry.S);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* skip error-related cols */
    c = skip_cols(c, err_cols);
    /* resolve t0 at pre-last column */
    const char *ptr;
    if (resolve_icgem_data_date_v2(skip_ws(c), entry.t0, ptr))
      ++error;
    /* resolve t1 at last column */
    if (!error)
      error += resolve_icgem_data_date_v2(skip_ws(ptr), entry.t1, ptr);
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
    auto res = std::from_chars(skip_ws(c + 4), line + len, entry.degree);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), line + len, entry.order);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), line + len, entry.C);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), line + len, entry.S);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* skip error-related cols */
    c = skip_cols(c, err_cols);
    /* resolve t0 at pre-last column */
    const char *ptr;
    if (resolve_icgem_data_date_v2(skip_ws(c), entry.t0, ptr))
      ++error;
    /* resolve t1 at last column */
    if (!error)
      error += resolve_icgem_data_date_v2(skip_ws(ptr), entry.t1, ptr);
    /* resolve period */
    res = std::from_chars(skip_ws(ptr), line + len, entry.period);
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
    auto res = std::from_chars(skip_ws(c + 4), line + len, entry.degree);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), line + len, entry.order);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), line + len, entry.C);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), line + len, entry.S);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    /* skip error-related cols */
    c = skip_cols(c, err_cols);
    /* resolve t0 at pre-last column */
    const char *ptr;
    if (resolve_icgem_data_date_v2(skip_ws(c), entry.t0, ptr))
      ++error;
    /* resolve t1 at last column */
    if (!error)
      error += resolve_icgem_data_date_v2(skip_ws(ptr), entry.t1, ptr);
    /* resolve period */
    res = std::from_chars(skip_ws(ptr), line + len, entry.period);
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
