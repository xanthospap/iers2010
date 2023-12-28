#include "aod1b.hpp"
#include <cstdio>

namespace {
inline const char *skipws(const char *line) noexcept {
  while (*line && *line == ' ')
    ++line;
  return line;
}
} /* unnamed namespace */

int dso::Aod1bIn::skip_header(std::ifstream &fin) const noexcept {
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Invalid stream for AOD1B file %s (traceback: %s)\n",
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

int dso::Aod1bIn::skip_lines(std::ifstream &fin, int num_lines) const noexcept {
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Invalid stream for AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  /* don't even get into trouble if lines are 0 */
  if (!num_lines) return 0;

  /* quick read requested number of lines */
  constexpr const int MAXL = 124;
  char line[MAXL];
  int linenr = 0;
  while (fin.getline(line, MAXL) && ++linenr < num_lines)

    if (!fin.good()) {
      if (fin.eof())
        return -1;
      return 1;
    }

  return 0;
}
