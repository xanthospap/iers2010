#include "aod1b.hpp"
#include "datetime/datetime_read.hpp"
#include <charconv>
#include <cstdio>
#include <cstring>

namespace {
inline const char *skipws(const char *line) noexcept {
  while (*line && *line == ' ')
    ++line;
  return line;
}
} /* unnamed namespace */

int dso::Aod1bIn::parse_data_block_header(
    const char *line, dso::Aod1bIn::Aod1bBlockHeader &rec) const noexcept {
  if (std::strncmp(line, "DATA SET", 8)) {
    return 1;
  }

  const char *str = line + 8;
  const char *last = line + std::strlen(line);
  int error = 0;

  auto res = std::from_chars(skipws(str), last, rec.mset_nr);
  error += (res.ec != std::errc{});
  str = res.ptr + 1;

  res = std::from_chars(skipws(str), last, rec.mnum_lines);
  error += (res.ec != std::errc{});
  str = res.ptr;

  if (std::strncmp(skipws(str), "COEFFICIENTS FOR ", 17))
    ++error;
  str = skipws(str + 17);
  rec.mepoch = dso::from_char<dso::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSS,
                              dso::nanoseconds>(str, &res.ptr);

  str = res.ptr + 1;
  if (std::strncmp(skipws(str), "OF TYPE", 7))
    ++error;
  str = skipws(str + 7);

  if (!std::strncmp(skipws(str), "atm", 3))
    rec.mtype = dso::AOD1BCoefficientType::ATM;
  else if (!std::strncmp(skipws(str), "ocn", 3))
    rec.mtype = dso::AOD1BCoefficientType::OCN;
  else if (!std::strncmp(skipws(str), "glo", 3))
    rec.mtype = dso::AOD1BCoefficientType::GLO;
  else if (!std::strncmp(skipws(str), "oba", 3))
    rec.mtype = dso::AOD1BCoefficientType::OBA;
  else
    ++error;

  if (error) {
    fprintf(
        stderr,
        "[ERROR] Failed parsing AOD1B block header line [%s] (tracback: %s)\n",
        line, __func__);
    return 1;
  }

  return 0;
}
