#include "aod1b.hpp"
#include <cstdio>

int dso::Aod1bIn::goto_next_block(
    std::ifstream &fin, dso::Aod1bIn::Aod1bBlockHeader &rec) const noexcept {
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Invalid stream for AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  constexpr const int MAXL = 124;
  const int MAXLS = num_data_sets();
  char line[MAXL];
  int error = 0, j = 0;
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

  if (!fin.good() || j >= MAXLS) {
    if (fin.eof())
      return -1;
    return 1;
  }

  return 0;
}

int dso::Aod1bIn::goto_next_block(
    std::ifstream &fin, dso::AOD1BCoefficientType type,
    dso::Aod1bIn::Aod1bBlockHeader &rec) const noexcept {
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Invalid stream for AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  constexpr const int MAXL = 124;
  const int MAXLS = num_data_sets();
  char line[MAXL];
  int error = 0, j = 0;
  while (fin.getline(line, MAXL) && (!error) && (++j < MAXLS)) {
    /* assuming next line was a data block header */
    if (parse_data_block_header(line, rec)) {
      ++error;
    } else {
    /* did we match the type ? */
    if (rec.mtype == type)
      return 0;
    /* skip lines untill next block */
    for (int l = 0; l < rec.mnum_lines; l++)
      fin.getline(line, MAXL);
    }
  }

  if (error) {
    fprintf(stderr,
            "[ERROR] Failed reading next header block for AOD1B file %s "
            "(traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  if (!fin.good() || j >= MAXLS) {
    if (fin.eof())
      return -1;
    return 1;
  }

  return 0;
}
