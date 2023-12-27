#include "aod1b.hpp"
#include "datetime/datetime_read.hpp"
#include <cctype>
#include <charconv>
#include <cstdio>
#include <cstring>
#include <fstream>

namespace {
const char *header_field(const char *line) noexcept {
  const char *c = line;
  while (*c && *c != ':')
    ++c;
  c = c + 1;
  while (*c && *c == ' ')
    ++c;
  return c;
}
const char *header_field(const char *line, int &sz) noexcept {
  /* start in reverse, find last char that is not ws */
  const char *c = line + std::strlen(line);
  while (*c && *c == ' ')
    --c;
  const char *last = (--c);
  /* find delimeter, i.e. ':' */
  while (*c && *c != ':')
    --c;
  /* first non ws character after delimeter */
  c = header_field(c);
  sz = (last - c) + 1;
  return c;
}
} /* unnamed namespace */

int dso::Aod1bIn::read_header() noexcept {
  std::ifstream fin(mfn.c_str());
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening AOD1B file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
    return 1;
  }

  constexpr const int MAXLS = 100;
  constexpr const int MAXL = 124;
  char line[MAXL];
  int error = 0, j = 0;
  while (fin.getline(line, MAXL) && (!error) && (++j < MAXLS)) {
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
                line, mfn.c_str(), __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "FILE FORMAT 0=BINARY 1=ASCII", 28)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 28), last, file_format()).ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, mfn.c_str(), __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "NUMBER OF HEADER RECORDS", 24)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 24), last, num_header_records())
              .ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, mfn.c_str(), __func__);
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
    } else if (!std::strncmp(line, "TIME EPOCH (GPS TIME)", 21)) {
      time_epoch() =
          dso::from_char<dso::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSS,
                         dso::nanoseconds>(header_field(line + 21));
    } else if (!std::strncmp(line, "TIME FIRST OBS(SEC PAST EPOCH)", 30)) {
      const char *last = line + std::strlen(line);
      double sec;
      auto res = std::from_chars(header_field(line + 30), last, sec);
      if (res.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, mfn.c_str(), __func__);
        ++error;
      }
      nanoseconds isec(static_cast<nanoseconds::underlying_type>(sec * 1e9));
      first_epoch() = time_epoch();
      first_epoch().add_seconds(isec);
      /* validate */
      while (!std::isdigit(*res.ptr))
        ++res.ptr;
      dso::Datetime<dso::nanoseconds> tmp =
          dso::from_char<dso::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSS,
                         dso::nanoseconds>(res.ptr);
      assert(tmp == first_epoch());
    } else if (!std::strncmp(line, "TIME LAST OBS(SEC PAST EPOCH)", 29)) {
      const char *last = line + std::strlen(line);
      double sec;
      auto res = std::from_chars(header_field(line + 29), last, sec);
      if (res.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, mfn.c_str(), __func__);
        ++error;
      }
      nanoseconds isec(static_cast<nanoseconds::underlying_type>(sec * 1e9));
      last_epoch() = time_epoch();
      last_epoch().add_seconds(isec);
      /* validate */
      while (!std::isdigit(*res.ptr))
        ++res.ptr;
      dso::Datetime<dso::nanoseconds> tmp =
          dso::from_char<dso::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSS,
                         dso::nanoseconds>(res.ptr);
      assert(tmp == last_epoch());
    } else if (!std::strncmp(line, "NUMBER OF DATA RECORDS", 22)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 22), last, num_data_records())
              .ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, mfn.c_str(), __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "PRODUCT START CREATE TIME(UTC)", 30)) {
      ;
    } else if (!std::strncmp(line, "PRODUCT END CREATE TIME(UTC)", 28)) {
      ;
    } else if (!std::strncmp(line, "FILESIZE (BYTES)", 16)) {
      ;
    } else if (!std::strncmp(line, "FILENAME", 8)) {
      ;
    } else if (!std::strncmp(line, "PROCESS LEVEL (1A OR 1B)", 24)) {
      int sz = 0;
      const char *r = header_field(line + 24, sz);
      std::memcpy(process_level(), r, sizeof(char) * sz);
    } else if (!std::strncmp(line, "PRESSURE TYPE (SP OR VI)", 24)) {
      int sz = 0;
      const char *r = header_field(line + 24, sz);
      std::memcpy(pressure_type(), r, sizeof(char) * sz);
    } else if (!std::strncmp(line, "MAXIMUM DEGREE", 14)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 14), last, max_degree()).ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, mfn.c_str(), __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "COEFFICIENTS ERRORS (YES/NO)", 28)) {
      int sz = 0;
      const char *r = header_field(line + 28, sz);
      if (sz == 2) {
        assert(!std::strncmp("NO", r, 2));
        coeff_errors() = 0;
      } else if (sz == 3) {
        assert(!std::strncmp("YES", r, 3));
        coeff_errors() = 1;
      } else {
        fprintf(
            stderr,
            "[ERROR] Unexpected \'COEFFICIENTS ERRORS\' header line in AOD1B "
            "file %s (traceback: %s)\n",
            mfn.c_str(), __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "COEFF. NORMALIZED (YES/NO)", 26)) {
      int sz = 0;
      const char *r = header_field(line + 26, sz);
      if (sz == 2) {
        assert(!std::strncmp("NO", r, 2));
        coeff_normalized() = 0;
      } else if (sz == 3) {
        assert(!std::strncmp("YES", r, 3));
        coeff_normalized() = 1;
      } else {
        fprintf(stderr,
                "[ERROR] Unexpected \'COEFF. NORMALIZED\' header line in AOD1B "
                "file %s (traceback: %s)\n",
                mfn.c_str(), __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "CONSTANT GM [M^3/S^2]", 21)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 21), last, GM()).ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, mfn.c_str(), __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "CONSTANT A [M]", 14)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 22), last, Re()).ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, mfn.c_str(), __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "CONSTANT FLAT [-]", 17)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 17), last, flat()).ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, mfn.c_str(), __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "CONSTANT OMEGA [RAD/S]", 22)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 17), last, omega()).ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, mfn.c_str(), __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "NUMBER OF DATA SETS", 19)) {
      const char *last = line + std::strlen(line);
      if (std::from_chars(header_field(line + 17), last, num_data_sets()).ec !=
          std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing line: [%s] from AOD1B file %s "
                "(traceback: %s)\n",
                line, mfn.c_str(), __func__);
        ++error;
      }
    } else if (!std::strncmp(line, "DATA FORMAT (N,M,C,S)", 21)) {
      ;
    } else if (!std::strncmp(line, "END OF HEADER", 13)) {
      break;
    } else {
      fprintf(stderr,
              "[ERROR] Unrecognized header line [%s] in AOD1B file %s "
              "(traceback: %s)\n",
              line, mfn.c_str(), __func__);
      ++error;
    }
  }

  /* hopefully no error occured with the stream! */
  if (!fin.good() || j >= MAXLS) {
    fprintf(
        stderr,
        "[ERROR] Failed reading header from AOD1B file %s (traceback: %s)\n",
        mfn.c_str(), __func__);
    return 1;
  }

  /* all done */
  return 0;
}
