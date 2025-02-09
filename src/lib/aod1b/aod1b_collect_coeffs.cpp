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

int dso::Aod1bIn::collect_coeffs(std::ifstream &fin, int max_degree,
                                 int max_order, int max_lines,
                                 dso::StokesCoeffs &cs,
                                 bool allow_resizing) const noexcept {
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Invalid stream for AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  if ((max_degree < max_order) || (max_degree > this->max_degree())) {
    fprintf(stderr,
            "[ERROR] Invalid max degree/order given for coefficient parsing; "
            "AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  if (!allow_resizing) {
    /* re-sizing not allowed, sizes must agree */
    if ((cs.max_degree() < max_degree) || (cs.max_order() < max_order)) {
      fprintf(
          stderr,
          "[ERROR]  Invalid max degree/order given for coefficient parsing; "
          "requested %dx%d but the StokesCoeffs instance can only hold %dx%d "
          "(traceback: %s)\n",
          max_degree, max_order, cs.max_degree(), cs.max_order(), __func__);
      return 1;
    }
  } else {
    /* re-sizing allowed; if needed set initial correct size of Stokes. note
     * that later we should again resize the instance to exactly match the 
     * size of coeffs (actually) collected.
     */
    cs.resize(max_degree, max_order);
  }

  /* set coeffs to zero */
  cs.clear();

  /* lines read/parsed */
  int lines_parsed = 0;

  constexpr const int MAXL = 124;
  char line[MAXL];
  int error = 0;
  /* actual degree/order collectedd from file */
  int max_degree_collected = 0;
  int max_order_collected = 0;

  /* watch the order in the while condition! */
  while ((lines_parsed < max_lines) && (!error) && fin.getline(line, MAXL)) {
    /* Older versions of gcc complain about uninitialized l, m values */ 
    int l = 0, m = 0;
    const char *stop = line + std::strlen(line);
    auto res = std::from_chars(skipws(line), stop, l);
    error += (res.ec != std::errc{});
    res = std::from_chars(skipws(res.ptr), stop, m);
    error += (res.ec != std::errc{});
    if ((l <= max_degree) && (m <= max_order) && (!error)) {
      res = std::from_chars(skipws(res.ptr), stop, cs.C(l, m));
      error += (res.ec != std::errc{});
      res = std::from_chars(skipws(res.ptr), stop, cs.S(l, m));
      error += (res.ec != std::errc{});
      max_degree_collected = std::max(max_degree_collected, l);
      max_order_collected = std::max(max_order_collected, m);
    }
    ++lines_parsed;
  }

  /* check for errors */
  if (error || !fin.good()) {
    fprintf(stderr,
            "[ERROR] Failed parsing coefficients for AOD1B file %s "
            "(traceback: %s)\n",
            mfn.c_str(), __func__);
    fprintf(stderr, "Parsing stopped at line [%s] i.e. %d/%d (traceback: %s)\n",
            line, lines_parsed, max_lines, __func__);
    return 1;
  }

  /* if resizing is allowed, then reset the Stokes coeffs to the right size */
  if (allow_resizing) {
    cs.cresize(max_degree_collected, max_order_collected);
  } else {
    /* resizing not allowed, we should have collected the exact number of
     * coefficients the user requested
     */
    if ((max_degree_collected != max_degree) ||
        (max_order_collected != max_order)) {
      fprintf(
          stderr,
          "[ERROR] Failed collecting the requested size of coefficients from "
          "AOD1B file! requested: %dx%d, collected: %dx%d (traceback: %s)\n",
          max_degree, max_order, max_degree_collected, max_order_collected,
          __func__);
      return -1;
    }
  }

  return 0;
}
