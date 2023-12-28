#include "aod1b.hpp"
#include <cctype>
#include <cstdio>
#include <cstring>
#include <fstream>

namespace {
inline const char *skipws(const char *line) noexcept {
  while (*line && *line == ' ') ++line;
  return line;
}
}/* unnamed namespace */

int dso::Aod1bIn::skip_header(std::ifstream &fin) const noexcept {
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Invalid stream for AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  fin.seekg(0, std::ios::beg);

  constexpr const int MAXL = 124;
  constexpr const int MAXLS = 100;
  char line[MAXL];
  int error = 0, j = 0;
  while (fin.getline(line, MAXL) && (!error) && (++j < MAXLS)) {
    if (!std::strncmp(line, "END OF HEADER", 13))
      return 0;
  }

  fprintf(
      stderr,
      "[ERROR] Failed finding END OF HEADER in AOD1B file %s (traceback: %s)\n",
      mfn.c_str(), __func__);
  return 1;
}

int dso::Aod1bIn::goto_next_block(std::ifstream &fin, dso::Aod1bIn::Aod1bBlockHeader &rec) const noexcept {
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Invalid stream for AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  constexpr const int MAXL = 124;
  const int MAXLS = num_data_sets();
  char line[MAXL];
  int error = 0, j=0;
  while (fin.getline(line, MAXL) && (!error) && (++j < MAXLS)) {
    if (parse_data_block_header(line, rec)) {
      ++error;
    }
  }
      
  if (error) {
      fprintf(stderr,
              "[ERROR] Failed reading next header block for AOD1B file %s "
              "(traceback: %s)\n",
              mfn.c_str(), __func__);
      return 1;
  }

  if (!fin.good() || j>=MAXLS) {
    if (fin.eof()) return -1;
    return 1;
  }

  return 0;
}

int dso::Aod1bIn::skip_lines(std::ifstream &fin, int num_lines) const noexcept {
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Invalid stream for AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  constexpr const int MAXL = 124;
  char line[MAXL];
  int linenr = 0;
  while (fin.getline(line, MAXL) && linenr++ < num_lines)

  if (!fin.good()) {
    if (fin.eof()) return -1;
    return 1;
  }

  return 0;
}

int dso::Aod1bIn::goto_next_block(
    std::ifstream &fin,
    dso::AOD1BCoefficientType type,
    dso::Aod1bIn::Aod1bBlockHeader &rec) const noexcept {
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Invalid stream for AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  constexpr const int MAXL = 124;
  const int MAXLS = num_data_sets();
  char line[MAXL];
  int error = 0, j=0;
  while (fin.getline(line, MAXL) && (!error) && (++j < MAXLS)) {
    if (parse_data_block_header(line, rec)) {
      ++error;
    }
    if (!error && rec.mtype == type)
      return 0;
    /* skip lines untill next block */
    for (int l=0; l<rec.mnum_lines; l++) fin.getline(line, MAXL);
  }
      
  if (error) {
      fprintf(stderr,
              "[ERROR] Failed reading next header block for AOD1B file %s "
              "(traceback: %s)\n",
              mfn.c_str(), __func__);
      return 1;
  }

  if (!fin.good() || j >= MAXLS) {
    if (fin.eof()) return -1;
    return 1;
  }

  return 0;
}
