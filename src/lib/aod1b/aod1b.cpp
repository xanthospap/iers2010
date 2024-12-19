#include "aod1b.hpp"
#include <cstdio>
#include <stdexcept>

dso::Aod1bIn::Aod1bIn(const char *fn) : mfn(fn) {
  std::memset(charArena, '\0', MaxArenaChars);
  if (read_header()) {
    fprintf(stderr,
            "[ERROR] Failed reading header off from AOD1B file %s (traceback: "
            "%s)\n",
            fn, __func__);
    throw std::runtime_error(
        "[ERROR] Failed reading header off from AOD1B file\n");
  }
}

dso::Aod1bIn::Aod1bIn(const dso::Aod1bIn &other) noexcept
    : mfn(other.mfn), mfile_type(other.mfile_type),
      mfile_format(other.mfile_format),
      mnum_header_records(other.mnum_header_records),
      mnum_data_records(other.mnum_data_records),
      mmax_degree(other.mmax_degree), mcoeff_errors(other.mcoeff_errors),
      mcoeff_normalized(other.mcoeff_normalized), mGM(other.mGM),
      mRe(other.mRe), mflat(other.mflat), momega(other.momega),
      mnum_data_sets(other.mnum_data_sets), mtime_epoch(other.mtime_epoch),
      mfirst_epoch(other.mfirst_epoch), mlast_epoch(other.mlast_epoch), mwave(other.mwave) {
  std::memcpy(charArena, other.charArena, sizeof(char) * 80);
}

dso::Aod1bIn::Aod1bIn(dso::Aod1bIn &&other) noexcept
    : mfn(std::move(other.mfn)), mfile_type(other.mfile_type),
      mfile_format(other.mfile_format),
      mnum_header_records(other.mnum_header_records),
      mnum_data_records(other.mnum_data_records),
      mmax_degree(other.mmax_degree), mcoeff_errors(other.mcoeff_errors),
      mcoeff_normalized(other.mcoeff_normalized), mGM(other.mGM),
      mRe(other.mRe), mflat(other.mflat), momega(other.momega),
      mnum_data_sets(other.mnum_data_sets), mtime_epoch(other.mtime_epoch),
      mfirst_epoch(other.mfirst_epoch), mlast_epoch(other.mlast_epoch), mwave(other.mwave) {
  std::memcpy(charArena, other.charArena, sizeof(char) * 80);
}

dso::Aod1bIn &dso::Aod1bIn::operator=(const dso::Aod1bIn &other) noexcept {
  if (this != &other) {
    mfn = other.mfn;
    mfile_type = other.mfile_type;
    mfile_format = other.mfile_format;
    mnum_header_records = other.mnum_header_records;
    mnum_data_records = other.mnum_data_records;
    mmax_degree = other.mmax_degree;
    mcoeff_errors = other.mcoeff_errors;
    mcoeff_normalized = other.mcoeff_normalized;
    mGM = other.mGM;
    mRe = other.mRe;
    mflat = other.mflat;
    momega = other.momega;
    mnum_data_sets = other.mnum_data_sets;
    mtime_epoch = other.mtime_epoch;
    mfirst_epoch = other.mfirst_epoch;
    mlast_epoch = other.mlast_epoch;
    mwave = other.mwave;
    std::memcpy(charArena, other.charArena, sizeof(char) * 80);
  }
  return *this;
}

dso::Aod1bIn &dso::Aod1bIn::operator=(dso::Aod1bIn &&other) noexcept {
  if (this != &other) {
    mfn = std::move(other.mfn);
    mfile_type = other.mfile_type;
    mfile_format = other.mfile_format;
    mnum_header_records = other.mnum_header_records;
    mnum_data_records = other.mnum_data_records;
    mmax_degree = other.mmax_degree;
    mcoeff_errors = other.mcoeff_errors;
    mcoeff_normalized = other.mcoeff_normalized;
    mGM = other.mGM;
    mRe = other.mRe;
    mflat = other.mflat;
    momega = other.momega;
    mnum_data_sets = other.mnum_data_sets;
    mtime_epoch = other.mtime_epoch;
    mfirst_epoch = other.mfirst_epoch;
    mlast_epoch = other.mlast_epoch;
    mwave = other.mwave;
    std::memcpy(charArena, other.charArena, sizeof(char) * 80);
  }
  return *this;
}
