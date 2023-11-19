#include "icgemio.hpp"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <charconv>

namespace {
constexpr int max_header_lines = 1000;
constexpr int bsz = 256;

/** Given a C-String, go to the first non-whitespace character. If the line is
 * empty or only contains whitespace characters, the returned (pointer to
 * char) is NULL
 * example: given line="  foo bar\0"
 *          return a pointer to 'f', aka the string "foo bar"
 * example: given "   \0"
 *          return '\0'
 */
const char *next_non_ws_char(const char *line) noexcept {
  const char *s = line;
  while (*s && *s == ' ')
    ++s;
  return s;
}

/** Given a C-String, go to the first whitespace character. */
const char *next_ws_char(const char *line) noexcept {
  const char *s = line;
  while (*s && *s != ' ')
    ++s;
  return s;
}

/** Given a null-terminated C-String, right trim whitespace characters */
const char *right_trim(char *str) noexcept {
  /* go to the rightmost character of the string */
  char *right = str + std::strlen(str) - 1;
  while (*right && right!=str && *right == ' ') {
    *right = '\0';
    --right;
  }
  return right;
}

/** Check if a given line corresponds to a given parameter name (or keyword)
 *
 * Example of ICGEM (header):
 * ---------------------------------------------------------------------------
 *           product_type gravity_field
 *              modelname EXAMPLE_MODEL
 * earth_gravity_constant  0.3986004415E+15
 *                 radius  0.6378136460E+07
 *             max_degree    370
 *                 errors formal
 * ---------------------------------------------------------------------------
 * So, to match a keyword, we should first skip any whitespace characters, and
 * then compare the first string encountered (up untill the next whitespace
 * character or EOL).
 *
 * @param[in] parameter Keyword to match (e.g. 'product_type', 'modelname',
 *            etc), as a null-terminated C-string.
 * @param[in] line Line to check (as in the header section)
 * @param[out] argument If the parameter (passed in) is indeed the keyword in
 *                 the passed in line, argument will point to the start of
 *                 the keyword argument (if any).
 *                 E.g., in line "max_degree    370", argument will point to
 *                 the character "3" (of 370).
 * @return True if the parameter is the line's keyword. Else return false.
 */
int match_header_keyword(const char *parameter, const char *line,
                         char const *&argument) noexcept {
  auto psz = std::strlen(parameter);

  /* go to the first word (skip whitespaces) */
  const char *str = next_non_ws_char(line);

  /* does the line start with the parameter_name string ? */
  if (!std::strncmp(str, parameter, psz)) {
    /* start of new word after the keyword in line */
    argument = next_non_ws_char(next_ws_char(str));
    return 1;
  } else {
    /* parameter_name not matched in begining of line; return NULL */
    return 0;
  }
}
} /* anonymous namespace */

/** @warning this function assumes that comment and header lines do not exceed
 *           256 characters; this is not guaranteed by the icgem format, but
 *           practically holds for all cases.
 */
int dso::Icgem::parse_header(bool quiet_skip) noexcept {

  /* open the ICGEM file */
  std::ifstream fin(_filename.c_str());
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening icgem file %s (traceback: %s)\n",
            _filename.c_str(), __func__);
    return 1;
  }

  char line[bsz];
  const char *arg;
  int counter = 0;
  int error = 0;

  /* first off, skip the part between the begining of the file and the start 
   * of the header. That is, find the line begining with the string:
   * 'begin_of_head'
   */
  fin.getline(line, bsz);
  while (std::strncmp(line, "begin_of_head", 13)) {
    if (!fin.getline(line, bsz) || (++counter > max_header_lines)) {
      fprintf(stderr,
              "[ERROR] Failed parsing header in icgem %s; too many header "
              "lines ... (traceback: %s)\n",
              _filename.c_str(), __func__);
      return 1;
    }
  }

  /* keep reading new lines untill we meet a line starting with 'end_of_head' */
  while (std::strncmp(line, "end_of_head", 11)) {

    /* check if this the 'product_type' field (mandatory keyword) */
    if (match_header_keyword("product_type", line, arg)) {
      if (!(*arg)) {
        fprintf(stderr,
                "[ERROR] Failed parsing value for parameter %s in icgem "
                "%s(traceback: %s)\n",
                "product_type", _filename.c_str(), __func__);
        error = 1;
      } else {
        right_trim(std::strcpy(_product_type, arg));
      }

      /* check if this the 'modelname' field (mandatory keyword) */
    } else if (match_header_keyword("modelname", line, arg)) {
      if (!(*arg)) {
        fprintf(stderr,
                "[ERROR] Failed parsing value for parameter %s in icgem "
                "%s(traceback: %s)\n",
                "modelname", _filename.c_str(), __func__);
        error = 1;
      } else {
        right_trim(std::strcpy(_modelname, arg));
      }

      /* check if this the 'earth_gravity_constant' field (mandatory keyword) */
    } else if (match_header_keyword("earth_gravity_constant", line, arg)) {
      if (!(*arg)) {
        fprintf(stderr,
                "[ERROR] Failed parsing value for parameter %s in icgem "
                "%s(traceback: %s)\n",
                "earth_gravity_constant", _filename.c_str(), __func__);
        error = 1;
      } else {
        if ((std::from_chars(arg, arg + std::strlen(arg),
                             _earth_gravity_constant))
                .ec != std::errc{}) {
          fprintf(stderr,
                  "[ERROR] Failed parsing value for parameter %s in icgem "
                  "%s(traceback: %s)\n",
                  "earth_gravity_constant", _filename.c_str(), __func__);
          error = 2;
        }
      }

      /* check if this the 'radius' field (mandatory keyword) */
    } else if (match_header_keyword("radius", line, arg)) {
      if (!(*arg)) {
        fprintf(stderr,
                "[ERROR] Failed parsing value for parameter %s in icgem "
                "%s(traceback: %s)\n",
                "radius", _filename.c_str(), __func__);
        error = 1;
      } else {
        if ((std::from_chars(arg, arg + std::strlen(arg), _radius)).ec !=
            std::errc{}) {
          fprintf(stderr,
                  "[ERROR] Failed parsing value for parameter %s in icgem "
                  "%s(traceback: %s)\n",
                  "radius", _filename.c_str(), __func__);
          error = 2;
        }
      }

      /* check if this the 'max_degree' field (mandatory keyword) */
    } else if (match_header_keyword("max_degree", line, arg)) {
      if (!(*arg)) {
        fprintf(stderr,
                "[ERROR] Failed parsing value for parameter %s in icgem "
                "%s(traceback: %s)\n",
                "max_degree", _filename.c_str(), __func__);
        error = 1;
      } else {
        if ((std::from_chars(arg, arg + std::strlen(arg), _max_degree)).ec !=
            std::errc{}) {
          fprintf(stderr,
                  "[ERROR] Failed parsing value for parameter %s in icgem "
                  "%s(traceback: %s)\n",
                  "max_degree", _filename.c_str(), __func__);
          error = 2;
        }
      }

      /* check if this the 'errors' field (mandatory keyword) */
    } else if (match_header_keyword("errors", line, arg)) {
      if (!(*arg)) {
        fprintf(stderr,
                "[ERROR] Failed parsing value for parameter %s in icgem "
                "%s(traceback: %s)\n",
                "errors", _filename.c_str(), __func__);
        error = 1;
      } else {
        right_trim(std::strcpy(_errors, arg));
      }

      /* check if this the 'tide_system' field (optional keyword) */
    } else if (match_header_keyword("tide_system", line, arg)) {
      if (!(*arg)) {
        fprintf(stderr,
                "[ERROR] Failed parsing value for parameter %s in icgem "
                "%s(traceback: %s)\n",
                "tide_system", _filename.c_str(), __func__);
        error = 1;
      } else {
        right_trim(std::strcpy(_tide_system, arg));
      }

      /* check if this the 'norm' field (optional keyword) */
    } else if (match_header_keyword("norm", line, arg)) {
      if (!(*arg)) {
        fprintf(stderr,
                "[ERROR] Failed parsing value for parameter %s in icgem "
                "%s(traceback: %s)\n",
                "norm", _filename.c_str(), __func__);
        error = 1;
      } else {
        right_trim(std::strcpy(_norm, arg));
      }

    } else {
      if (!quiet_skip) {
        fprintf(
            stderr,
            "[WARNING] Skipping header line %s of ICGEM file %s (traceback: "
            "%s)\n",
            line, _filename.c_str(), __func__);
      }
    }

    fin.getline(line, bsz);
    if (error || (++counter > max_header_lines)) {
      if (error)
        return error;
      fprintf(stderr,
              "[ERROR] Failed parsing header in icgem %s; too many header "
              "lines ... (traceback: %s)\n",
              _filename.c_str(), __func__);
      return 1;
    }

  } /* keep on parsing (header) lines ... */

  data_section_pos = fin.tellg();
  return 0;
}
