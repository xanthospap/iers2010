#include "tropo.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>

constexpr const int MAX_LINE = 512 * 2;

const char *HeaderStr =
    "%  lat    lon   p:a0    A1   B1   A2   B2  T:a0   A1    B1   A2   B2  "
    "Q:a0   A1    B1    A2    B2    dT:a0   A1   B1   A2   B2  undu      Hs  "
    "  a_h:a0     A1       B1        A2      B2     a_w:a0     A1       B1   "
    "    A2       B2    lambda:a0  A1     B1      A2     B2    Tm:a0   A1   "
    "B1   A2   B2    Gn_h:a0   A1    B1     A2     B2   Ge_h:a0   A1    B1   "
    "  A2     B2   Gn_w:a0   A1    B1     A2     B2   Ge_w:a0   A1    B1     "
    "A2     B2";

/// @brief Resolve a GPT3 grid file longitude resolution using the first two
///        lines
/// @param[in] fin An open GPT3 grid file (at whatever position)
/// @param[out] resolution Longitude resolution in [deg]. Either 1 or 5
/// @return Anything other than 0 denotes an error
int inspect_grid_file_resolution(std::ifstream &fin,
                                 double &resolution) noexcept {
  if (!fin.is_open() || !fin.good()) {
    return 1;
  }

  fin.seekg(0);

  int error;

  // get the first line, should be the header
  char line[MAX_LINE];
  if (fin.getline(line, MAX_LINE)) {
    unsigned sz = std::strlen(HeaderStr);
    error = std::strncmp(line, HeaderStr, sz);
  } else {
    error = 1;
  }
  if (error)
    return error;

  double lat1, lat2, lon1, lon2;
  char *start, *end;

  // get first data line
  fin.getline(line, MAX_LINE);
  start = line;
  lat1 = std::strtod(start, &end);
  if (lat1 == 0e0 || end == start)
    return 2;
  start = end;
  lon1 = std::strtod(start, &end);
  if (lon1 == 0e0 || end == start)
    return 2;

  // get second data line
  fin.getline(line, MAX_LINE);
  start = line;
  lat2 = std::strtod(start, &end);
  if (lat2 == 0e0 || end == start)
    return 2;
  start = end;
  lon2 = std::strtod(start, &end);
  if (lon2 == 0e0 || end == start)
    return 2;

  // resolve
  if (lat1 == lat2) {
    const int dlon = lon2 * 10 - lon1 * 10;
    if (dlon == 10)
      resolution = 1e0;
    else if (dlon == 50)
      resolution = 5e0;
    else
      error = 3;
  } else {
    error = 3;
  }

  return error;
}

int dso::Gpt3Grid::parse_grid(const char *gridfn) noexcept {
  // open the grid file
  std::ifstream grd(gridfn);
  if (!grd) {
    fprintf(stderr,
            "[ERROR] Failed to open gpt3 grid file %s (traceback: %s)\n",
            gridfn, __func__);
    return 1;
  }

  // check the resolution
  double resolution;
  if (int error; (error=inspect_grid_file_resolution(grd, resolution)) ) {
    fprintf(stderr,
            "[ERROR] Failed to find valid resolution for gpt3 grid file %s, error nr:%d "
            "(traceback: %s)\n",
            gridfn,error, __func__);
    return 1;
  }

  // allocate memmory (if needed)
  if (allocate(resolution))
    return 1;

  char line[MAX_LINE];

  // go to top of file and read/skip header. header is skiped because it is
  // already validated in the inspect_grid_file_resolution function
  grd.seekg(0);
  grd.getline(line, MAX_LINE);

  // parse the file ...
  int row = 0, error = 0;
  while (grd.getline(line, MAX_LINE) && !error) {
    char *start = line, *end;

    [[maybe_unused]] double dummy = std::strtod(start, &end);

    start = end;
    dummy = std::strtod(start, &end);

    start = end;
    double *data = this->p_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end);
      error += (start == end);
      start = end;
    }

    data = this->t_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end);
      error += (start == end);
      start = end;
    }

    data = this->q_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end) * 1e-3;
      error += (start == end);
      start = end;
    }

    data = this->dt_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end) * 1e-3;
      error += (start == end);
      start = end;
    }

    // data = this->u_data_row();
    *(u_grid_row() + row) = std::strtod(start, &end);
    error += (start == end);
    start = end;

    //data = this->hs_data_row();
    *(hs_grid_row() + row) = std::strtod(start, &end);
    error += (start == end);
    start = end;

    data = this->ah_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end) * 1e-3;
      error += (start == end);
      start = end;
    }

    data = this->aw_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end) * 1e-3;
      error += (start == end);
      start = end;
    }

    data = this->la_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end);
      error += (start == end);
      start = end;
    }

    data = this->tm_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end);
      error += (start == end);
      start = end;
    }

    data = this->gn_h_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end) * 1e-5;
      error += (start == end);
      start = end;
    }

    data = this->ge_h_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end) * 1e-5;
      error += (start == end);
      start = end;
    }

    data = this->gn_w_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end) * 1e-5;
      error += (start == end);
      start = end;
    }

    data = this->ge_w_grid_row(row);
    for (int i = 0; i < 5; i++) {
      data[i] = std::strtod(start, &end) * 1e-5;
      error += (start == end);
      start = end;
    }

    ++row;
  }

  if (error) {
    fprintf(stderr,
            "[ERROR] Failed parsing gpt3 grid file %s; line nr %d (traceback: "
            "%s)\n",
            gridfn, row + 2, __func__);
    return 1;
  }

  if (row != num_rows) {
    fprintf(stderr,
            "[ERROR] Failed parsing gpt3 grid file %s; read %d/%d data lines "
            "(traceback: %s)\n",
            gridfn, row, num_rows, __func__);
    return 1;
  }

  return 0;
}
