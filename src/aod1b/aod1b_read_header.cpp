#include "aod1b.hpp"
#include "datetime/datetime_read.hpp"
#include <cstdio>
#include <charconv>
#include <fstream>
#include <cstring>
#include  <cctype>

namespace {
const char *header_field(const char *line, int &sz) noexcept {
  /* start in reverse, find last char that is not ws */
  const char *c = line + std::strlen(line);
  while (*c && *c == ' ') --c;
  const char *last = c;
  while (*c && *c != ':') --c;
  while (*c && *c == ' ') ++c;
  sz = last - c;
  return c;
}
const char *header_field(const char *line) noexcept {
  const char *c = line;
  while (*c && *c != ':') ++c;
  while (*c && *c == ' ') ++c;
  return c;
}
}/* unnamed namespace */

int dso::Aod1bIn::read_header(const char *fn) noexcept {
  std::ifstream fin(fn);
  if (!ifstream.open()) {
    fprintf(stderr, "[ERROR] Failed opening AOD1B file %s (traceback: %s)\n", fn, __func__);
    return 1;
  }

  constexpr const int MAXL = 124;
  char line[MAXL];
  int error = 0;
  while (fin.getline(line, MAXL) && (!error)) {
    if (!std::strncmp(line, "PRODUCER AGENCY", 15)) {
      int sz = 0;
      const char *r = header_field(line + 15, sz);
      std::memcpy(agency(), r, sizeof(char) * sz);
    } else if (!std::strncmp(line, "PRODUCER INSTITUTION", 20)) {
      /* do nothing */
      ;
    } else if (!std::strncmp(line, "FILE TYPE ipAOD1BF", 18)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 18), last, file_type()).ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, fn, __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "FILE FORMAT 0=BINARY 1=ASCII", 28)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 28), last, file_format()).ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, fn, __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "NUMBER OF HEADER RECORDS", 24)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 24), last, num_header_records()).ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, fn, __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "SOFTWARE VERSION", 16)) {
      ;
    } else if (!std::strncmp(line, "SOFTWARE LINK TIME", 18)) {
      ;
    } else if (!std::strncmp(line, "REFERENCE DOCUMENTATION", 23)) {
      ;
    } else if (!std::strncmp(line, "SATELLITE NAME", 14)) {
      int sz = 0;
      const char *r = header_field(line + 14, sz);
      std::memcpy(satellite(), r, sizeof(char) * sz);
    } else if (!std::strncmp(line, "SENSOR NAME", 11)) {
      int sz = 0;
      const char *r = header_field(line + 11, sz);
      std::memcpy(sensor(), r, sizeof(char) * sz);
    } else if (!std::strncmp(line, "TIME EPOCH (GPS TIME)", 23)) {
      time_epoch() =
          dso::from_char<ods::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSS,
                         dso::nanoseconds>(header_field(line + 23));
    } else if (!std::strncmp(line, "TIME FIRST OBS(SEC PAST EPOCH)", 30)) {
      const char *last = line + std::strlen(line);
      double sec;
      auto res = std::from_chars(header_field(line + 30), last, sec);
      if (res.ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, fn, __func__);
        ++error;
      }
      nanoseconds isec(static_cast<nanoseconds::underlying_type>(sec * 1e9));
      first_epoch() = time_epoch();
      first_epoch().add_seconds(isec);
      /* validate */
      while (!std::isdigit(res.ptr)) ++res.ptr;
      dso::Datetime<dso::nanoseconds> tmp =
          dso::from_char<dso::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSS,
                         dso::nanoseconds>(res.ptr);
      assert(tmp == first_epoch());
    } else if (!std::strncmp(line, "TIME LAST OBS(SEC PAST EPOCH)", 29)) {
      const char *last = line + std::strlen(line);
      double sec;
      auto res = std::from_chars(header_field(line + 29), last, sec);
      if (res.ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, fn, __func__);
        ++error;
      }
      nanoseconds isec(static_cast<nanoseconds::underlying_type>(sec * 1e9));
      last_epoch() = time_epoch();
      last_epoch().add_seconds(isec);
      /* validate */
      while (!std::isdigit(res.ptr)) ++res.ptr;
      dso::Datetime<dso::nanoseconds> tmp =
          dso::from_char<dso::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSS,
                         dso::nanoseconds>(res.ptr);
      assert(tmp == last_epoch());
    } else if (!std::strncmp(line, "NUMBER OF DATA RECORDS", 22)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 22), last, num_data_records()).ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, fn, __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "PRODUCT START CREATE TIME(UTC)", 30)) {
      ;
    } else if (!std::strncmp(line, "PRODUCT END CREATE TIME(UTC)", 28)) {
      ;
    } else if (!std::strncmp(line, "FILESIZE (BYTES)", 16)) {
      ;
    } else if (!std::strncmp(line, "FILESIZE (BYTES)", 16)) {
      ;
  }
}
