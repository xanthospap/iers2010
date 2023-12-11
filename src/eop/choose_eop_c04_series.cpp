#include "eop.hpp"
#include <charconv>
#include <cstdio>
#include <cstring>
#include <fstream>

namespace {
constexpr const std::size_t MAX_LINE_CHARS = 256;
const char *h1c0420 = "# EARTH ORIENTATION PARAMETER (EOP) PRODUCT CENTER "
                      "CENTER (PARIS OBSERVATORY)";
const char *h2c0420 =
    "# EOP (IERS) 20 C04 TIME SERIES  consistent with ITRF 2020 - sampled";

const char *h1c0414 = "EARTH ORIENTATION PARAMETER (EOP) PRODUCT CENTER CENTER "
                      "(PARIS OBSERVATORY)";
const char *h2c0414 =
    "INTERNATIONAL EARTH ROTATION AND REFERENCE SYSTEMS SERVICE";
const char *h3c0414 = "EOP (IERS) 14 C04 TIME SERIES";

const char *remove_leading_ws(const char *line) noexcept {
  while (*line && *line == ' ')
    ++line;
  return line;
}
} /* unnamed namespace */

/*
 * Criteria C04/20:
 * [# EARTH ORIENTATION PARAMETER (EOP) PRODUCT CENTER CENTER (PARIS
 * OBSERVATORY)[...]
 * [# EOP (IERS) 20 C04 TIME SERIES consistent with ITRF 2020 - sampled[...]
 *
 * Criteria C04/14:
 * [               EARTH ORIENTATION PARAMETER (EOP) PRODUCT CENTER CENTER
 * (PARIS OBSERVATORY)[...] [                      INTERNATIONAL EARTH ROTATION
 * AND REFERENCE SYSTEMS SERVICE[...] [                                    EOP
 * (IERS) 14 C04 TIME SERIES[...]
 */

int dso::details::choose_c04_series(
    const char *c04fn, dso::details::IersEopFormat &type) noexcept {
  /* open file */
  std::ifstream fin(c04fn);
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed opening EOP (C04) file %s (traceback: "
            "%s)\n",
            c04fn, __func__);
    return 1;
  }

  char line[MAX_LINE_CHARS];
  fin.getline(line, MAX_LINE_CHARS);

  /* try for C04/20 */
  if (!std::strncmp(line, h1c0420, std::strlen(h1c0420))) {
    fin.getline(line, MAX_LINE_CHARS);
    if (!std::strncmp(line, h2c0420, std::strlen(h2c0420))) {
      type = dso::details::IersEopFormat::C0420;
      return 0;
    } else {
      fprintf(stderr,
              "[ERROR] Failed matching second header line in C04/20 file "
              "%s\n[ERROR] Line was [%s] (traceback: %s)\n",
              c04fn, line, __func__);
      return 1;
    }
  }

  /* try for C04/14 */
  if (!std::strncmp(remove_leading_ws(line), h1c0414, std::strlen(h1c0414))) {
    fin.getline(line, MAX_LINE_CHARS);
    if (std::strncmp(remove_leading_ws(line), h2c0414, std::strlen(h2c0414))) {
      fprintf(stderr,
              "[ERROR] Failed matching second header line in C04/14 file "
              "%s\n[ERROR] Line was [%s] (traceback: %s)\n",
              c04fn, line, __func__);
      return 1;
    }
    fin.getline(line, MAX_LINE_CHARS);
    if (std::strncmp(remove_leading_ws(line), h3c0414, std::strlen(h3c0414))) {
      fprintf(stderr,
              "[ERROR] Failed matching third header line in C04/20 file "
              "%s\n[ERROR] Line was [%s] (traceback: %s)\n",
              c04fn, line, __func__);
      return 1;
    } else {
      type = dso::details::IersEopFormat::C0414;
      return 0;
    }
  }

  /* no match for first line */
  fprintf(stderr,
          "[ERROR] Failed to match first (header) line in EOP file %s "
          "(traceback: %s)\n",
          c04fn, __func__);
  return 1;
}
