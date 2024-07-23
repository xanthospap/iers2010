#include "pole_tide.hpp"
#include <charconv>
#include <cstdio>
#include <cstring>
#include <fstream>

namespace {
inline bool cmp_nwsc(const char *line) noexcept {
  const char *header = "nmAnm(Real)Bnm(Real)Anm(Imaginary)Bnm(Imaginary)";
  const int hsz = std::strlen(header);
  const char *c = line;
  int idx = 0;
  while (*c) {
    if (*c != ' ') {
      if (*c == header[idx]) {
        ++idx;
      } else {
        return 1;
      }
    }
    ++c;
  }
  return (idx == hsz) ? 0 : 1;
}

inline const char *skipws(const char *line) noexcept {
  while (*line && *line == ' ')
    ++line;
  return line;
}
} /* unnamed namespace */

int dso::pole_tide_details::parse_desai02_coeffs(
    const char *fn, int max_degree, int max_order,
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> &A_real,
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> &A_imag,
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> &B_real,
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise>
        &B_imag) noexcept {

  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed to open file %s for reading Ocean Pole Tide "
            "coeffs (traceback: %s)\n",
            fn, __func__);
    return 1;
  }

  if ((max_degree < max_order) ||
      (max_degree > dso::pole_tide_details::MAX_DEGREE_DESAI_2002) ||
      (max_order > dso::pole_tide_details::MAX_ORDER_DESAI_2002)) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for parsing Ocean Pole Tide %s "
            "(traceback: %s)\n",
            fn, __func__);
    return 1;
  }

  /* resize (if needed) */
  A_real.resize(max_degree, max_order);
  A_imag.resize(max_degree, max_order);
  B_real.resize(max_degree, max_order);
  B_imag.resize(max_degree, max_order);

  constexpr const int SZ = 256;
  char line[SZ];

  /* read header and validate */
  fin.getline(line, SZ);
  if (cmp_nwsc(line)) {
    fprintf(stderr,
            "[ERROR] Failed validating header in Ocean Pole Tide coeffs file "
            "%s (traceback: %s)\n",
            fn, __func__);
    return 1;
  }

  int error = 0;
  int l, m;
  double data[4];
  while (fin.getline(line, SZ) && !error) {
    auto sz = std::strlen(line);
    const char *from = line;
    auto res = std::from_chars(skipws(from), line + sz, l);
    error += (res.ec != std::errc{});
    from = res.ptr;
    res = std::from_chars(skipws(from), line + sz, m);
    error += (res.ec != std::errc{});
    from = res.ptr;
    for (int i = 0; i < 4; i++) {
      res = std::from_chars(skipws(from), line + sz, data[i]);
      error += (res.ec != std::errc{});
      from = res.ptr;
    }
    if (l <= max_degree && m <= max_order) {
      A_real(l, m) = data[0];
      B_real(l, m) = data[1];
      A_imag(l, m) = data[2];
      B_imag(l, m) = data[3];
    }
  }

  if (error) {
    fprintf(stderr,
            "[ERROR] Failed resolving data in Ocean Pole Tide coeffs file %s "
            "(traceback: %s)\n",
            fn, __func__);
    fprintf(stderr, "[ERROR] Error line \"%s\" (traceback: %s)\n", line,
            __func__);
    return 1;
  }

  if (!fin.good() && !fin.eof()) {
    fprintf(stderr,
            "[ERROR] Failed resolving data in Ocean Pole Tide coeffs file %s "
            "(traceback: %s)\n",
            fn, __func__);
    fprintf(stderr, "[ERROR] Error line \"%s\" (traceback: %s)\n", line,
            __func__);
    return 5;
  }

  return 0;
}
