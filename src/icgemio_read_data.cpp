#ifdef FOOBAR

#include "datetime/dtcalendar.hpp"
#include "icgemio.hpp"
#include <cassert>
#include <charconv>

namespace {

const char *skip_ws(const char *line) noexcept {
  while (*line && *line == ' ') ++line;
  return line;
}

int num_columns(const char *line) noexcept {
  int num_cols = (*line != ' ' && *line);
  while (*(line + 1)) {
    if (*line == ' ' && *(line + 1) != ' ')
      ++num_cols;
    ++line;
  }
  return num_cols;
}

const char *goto_last_column(const char *line, int max_size = -1) noexcept {
  if (max_size < 0) max_size = std::strlen(line);
  if (!max_size) return nullptr;

  const char *s = line + max_size - 1;
  /* first right trim */
  while (s != line && *s == ' ')
    --s;
  if (s == line) return nullptr;
  /* now start looking for a whitespace character, from the top-right point */
  while (*s != ' ' && s >= line)
    --s;
  return s + 1;
}

const char *goto_last_two_columns(const char *line) noexcept {
  const char *s = goto_last_column(line);
  if (!s) return nullptr;
  int max_size = s - line;
  return goto_last_column(line, max_size);
}

const char *skip_cols(const char *line, int num_cols) noexcept {
  const char *str = skip_ws(line);
  int cols_skipped = 0;
  while (*str && *(str+1) && (cols_skipped<num_cols)) {
    if (*str != ' ' && *(str+1) == ' ') ++cols_skipped;
    ++str;
  }
  return skip_ws(str);
}

int resolve_data_line_v1(const char *line, DataEntry &entry) noexcept {
  int error = 0;
  int len = std::strlen(line);
  const char *c = line;

  /* Resolve Static spherical harmonic coefficients 
   * We are only interested in the parameters: l,m,C,S which are always the 
   * first four parameters (in this sequence).
   */
  if (!std::strncmp(line, "gfc ", 4)) {
    /* resolve degree */
    auto res = std::from_chars(skip_ws(c+4), line+len, entry.degree);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), line+len, entry.order);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), line+len, entry.C);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), line+len, entry.S);
    if (res.ec != std::errc{}) ++error;
  } else if (!std::strncmp(line, "gfct", 4)) {
    /* Time variable spherical harmonic coefficients, according to the 
     * respective formulae.
     *
     * degree, order, Clm, Slm, are always the first four columns; 
     * time (yyyymmdd) is always the last column.
     */
    /* resolve degree */
    auto res = std::from_chars(skip_ws(c+4), line+len, entry.degree);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), line+len, entry.order);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), line+len, entry.C);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), line+len, entry.S);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve t at last column */
    res =
        std::from_chars(skip_ws(goto_last_column(c)), line + len, entry.t0);
    if (res.ec != std::errc{}) ++error;
  } else if (!std::strncmp(line, "trnd", 4)) {
    /* Time variable spherical harmonic coefficients, according to the 
     * respective formulae.
     *
     * degree, order, Clm, Slm, are always the first four columns; 
     */
    /* resolve degree */
    auto res = std::from_chars(skip_ws(c+4), line+len, entry.degree);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), line+len, entry.order);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), line+len, entry.C);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), line+len, entry.S);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
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
    auto res = std::from_chars(skip_ws(c+4), line+len, entry.degree);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), line+len, entry.order);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), line+len, entry.C);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), line+len, entry.S);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve t at last column */
    res = std::from_chars(skip_ws(goto_last_column(c)), line + len,
                          entry.period);
    if (res.ec != std::errc{}) ++error;
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
    auto res = std::from_chars(skip_ws(c+4), line+len, entry.degree);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve order */
    res = std::from_chars(skip_ws(c), line+len, entry.order);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve C(degree, order) */
    res = std::from_chars(skip_ws(c), line+len, entry.C);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve S(degree, order) */
    res = std::from_chars(skip_ws(c), line+len, entry.S);
    if (res.ec != std::errc{}) ++error;
    c = res.ptr;
    /* resolve t at last column */
    res = std::from_chars(skip_ws(goto_last_column(c)), line + len,
                          entry.period);
    if (res.ec != std::errc{}) ++error;
  } else {
    fprintf(stderr,
            "[ERROR] Unknown keyword encountered in the ICGEM file! line [%s], "
            "(traceback: %s)\n",
            line, __func__);
    ++error;
  }
  return error;
}
} /* anonymous namespace */

/** @warning coeffs should have already been initialized and allocated with
 *           enough memmory to hold the (to-be-) parsed coefficients.
 */
int dso::Icgem::parse_data(int l, int m,
                           const dso::TwoPartDate &t,
                           dso::StokesCoeffs *coeffs) noexcept {

  /* clear out coeffs (i.e. set to zero) */
  coeffs->clear();

  /* t in fractional years (needed for TVG terms) */
  const double tyears = t.as_fractional_years();

#endif
