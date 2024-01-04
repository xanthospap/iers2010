#include "aod1b.hpp"
#include <cstdio>
#include <charconv>

namespace {
inline const char *skipws(const char *line) noexcept {
  while (*line && *line == ' ')
    ++line;
  return line;
}
} /* unnamed namespace */

int dso::Aod1bIn::get_tidal_wave(std::ifstream &fin, int max_degree,
                                 int max_order,
                                 dso::StokesCoeffs &cos_cs, dso::StokesCoeffs &sin_cs) const noexcept {
  if ((max_degree < max_order) || (max_degree > this->max_degree())) {
    fprintf(stderr,
            "[ERROR] Invalid max degree/order given for coefficient parsing; "
            "AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Invalid stream for AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  /* set coeffs to zero */
  cos_cs.clear();
  sin_cs.clear();

  constexpr const int MAXL = 124;
  char line[MAXL];
  int error = 0;
  int max_lines;
  int lines_parsed = 0;
  
  /* this should be the cos coeffs */
  fin.getline(line, MAXL);
  const char *end = line + std::strlen(line);
  if (std::strncmp(skipws(line), "DATA SET 01:", 12))
    ++error;
  if (std::from_chars(skipws(line), end, max_lines).ec != std::errc{})
    ++error;
  end += 20;
  if (std::strncmp(skipws(end), "COEFFICIENTS OF TYPE cos", 24))
    ++error;
  if (error) {
    fprintf(stderr, "[ERROR] failed resolving (tidal) header from line %s from AOD1B file %s (trceback: %s)\n", line, mfn.c_str(), __func__);
    return 1;
  }

  /* watch the order in the while condition! */
  while ((lines_parsed < max_lines) && (!error) && fin.getline(line, MAXL)) {
    int l, m;
    const char *stop = line + std::strlen(line);
    auto res = std::from_chars(skipws(line), stop, l);
    error += (res.ec != std::errc{});
    res = std::from_chars(skipws(res.ptr), stop, m);
    error += (res.ec != std::errc{});
    if ((l <= max_degree) && (m <= max_order) && (!error)) {
      res = std::from_chars(skipws(res.ptr), stop, cos_cs.C(l, m));
      error += (res.ec != std::errc{});
      res = std::from_chars(skipws(res.ptr), stop, cos_cs.S(l, m));
      error += (res.ec != std::errc{});
    }
    ++lines_parsed;
  }

  /* check for errors */
  if (error || !fin.good()) {
    fprintf(stderr,
            "[ERROR] Failed parsing coefficients for AOD1B file %s "
            "(traceback: %s)\n",
            mfn.c_str(), __func__);
    fprintf(stderr, "Parsing stopped at line [%s]\n", line);
    return 1;
  }
  
  /* this should be the sin coeffs */
  lines_parsed = 0;
  fin.getline(line, MAXL);
  const char *end = line + std::strlen(line);
  if (std::strncmp(skipws(line), "DATA SET 01:", 12))
    ++error;
  if (std::from_chars(skipws(line), end, max_lines).ec != std::errc{})
    ++error;
  end += 20;
  if (std::strncmp(skipws(end), "COEFFICIENTS OF TYPE sin", 24))
    ++error;
  if (error) {
    fprintf(stderr, "[ERROR] failed resolving (tidal) header from line %s from AOD1B file %s (trceback: %s)\n", line, mfn.c_str(), __func__);
    return 1;
  }

  /* watch the order in the while condition! */
  while ((lines_parsed < max_lines) && (!error) && fin.getline(line, MAXL)) {
    int l, m;
    const char *stop = line + std::strlen(line);
    auto res = std::from_chars(skipws(line), stop, l);
    error += (res.ec != std::errc{});
    res = std::from_chars(skipws(res.ptr), stop, m);
    error += (res.ec != std::errc{});
    if ((l <= max_degree) && (m <= max_order) && (!error)) {
      res = std::from_chars(skipws(res.ptr), stop, sin_cs.C(l, m));
      error += (res.ec != std::errc{});
      res = std::from_chars(skipws(res.ptr), stop, sin_cs.S(l, m));
      error += (res.ec != std::errc{});
    }
    ++lines_parsed;
  }

  /* check for errors */
  if (error || !fin.good()) {
    fprintf(stderr,
            "[ERROR] Failed parsing coefficients for AOD1B file %s "
            "(traceback: %s)\n",
            mfn.c_str(), __func__);
    fprintf(stderr, "Parsing stopped at line [%s]\n", line);
    return 1;
  }

  return 0;
}
