#include "costg_utils.hpp"
#include "eop.hpp"
#include "geodesy/units.hpp"
#include <charconv>
#include <cstdio>
#include <fstream>
#include <stdexcept>

const char *costg::skipws(const char *line) noexcept {
  while (line && *line == ' ')
    ++line;
  return line;
}

std::vector<costg::BmAcceleration> costg::parse_acceleration(const char *fn) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed opening file %s\n", fn);
    throw std::runtime_error("Failed opening input file " + std::string(fn) +
                             "\n");
  }

  char line[256];
  // for (int i = 0; i < 6; i++)
  //   fin.getline(line, 256);
  std::vector<costg::BmAcceleration> vec;

  while (fin.getline(line, 256)) {
    int error = 0;
    double td[4];
    int sz = std::strlen(line);
    const char *str = line;
    for (int i = 0; i < 4; i++) {
      auto res = std::from_chars(skipws(str), line + sz, td[i]);
      if (res.ec != std::errc{})
        ++error;
      str = res.ptr;
    }
    if (!error) {
      int imjd = (int)td[0];
      double fsec = (td[0] - imjd) * 86400e0;
      const dso::MjdEpoch t(imjd, dso::FractionalSeconds{fsec});
      vec.emplace_back(t, td[1], td[2], td[3]);
    }
  }

  if (!fin.eof()) {
    fprintf(stderr, "ERROR Failed parsing input file %s\n", fn);
    throw std::runtime_error("Failed parsing input file " + std::string(fn) +
                             "\n");
  }

  return vec;
}

std::vector<costg::BmOrbit> costg::parse_orbit(const char *fn) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed opening file %s\n", fn);
    throw std::runtime_error("Failed opening input file " + std::string(fn) +
                             "\n");
  }

  char line[512];
  // for (int i = 0; i < 6; i++)
  //   fin.getline(line, 512);
  std::vector<costg::BmOrbit> vec;

  while (fin.getline(line, 512)) {
    int error = 0;
    double td[10];
    int sz = std::strlen(line);
    const char *str = line;
    for (int i = 0; i < 10; i++) {
      auto res = std::from_chars(skipws(str), line + sz, td[i]);
      if (res.ec != std::errc{})
        ++error;
      str = res.ptr;
    }
    if (!error) {
      int imjd = (int)td[0];
      double fsec = (td[0] - imjd) * 86400e0;
      const dso::MjdEpoch t(imjd, dso::FractionalSeconds{fsec});
      vec.emplace_back(t, &td[1]);
    }
  }

  if (!fin.eof()) {
    fprintf(stderr, "ERROR Failed parsing input file %s\n", fn);
    throw std::runtime_error("Failed parsing input file " + std::string(fn) +
                             "\n");
  }

  return vec;
}

std::vector<costg::BmRotaryMatrix> costg::parse_rotary(const char *fn) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed opening file %s\n", fn);
    throw std::runtime_error("Failed opening input file " + std::string(fn) +
                             "\n");
  }

  char line[1024];
  // for (int i = 0; i < 6; i++)
  //   fin.getline(line, 1024);
  std::vector<costg::BmRotaryMatrix> vec;

  while (fin.getline(line, 1024)) {
    int error = 0;
    double td[10];
    int sz = std::strlen(line);
    const char *str = line;
    for (int i = 0; i < 10; i++) {
      auto res = std::from_chars(skipws(str), line + sz, td[i]);
      if (res.ec != std::errc{})
        ++error;
      str = res.ptr;
    }
    if (!error) {
      int imjd = (int)td[0];
      double fsec = (td[0] - imjd) * 86400e0;
      const dso::MjdEpoch t(imjd, dso::FractionalSeconds{fsec});
      vec.emplace_back(t, &td[1]);
    }
  }

  if (!fin.eof()) {
    fprintf(stderr, "ERROR Failed parsing input file %s\n", fn);
    throw std::runtime_error("Failed parsing input file " + std::string(fn) +
                             "\n");
  }

  return vec;
}

std::vector<costg::BmEops> costg::parse_eops(const char *fn) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed opening file %s\n", fn);
    throw std::runtime_error("Failed opening input file " + std::string(fn) +
                             "\n");
  }

  char line[1024];
  // for (int i = 0; i < 6; i++)
  //   fin.getline(line, 1024);
  std::vector<costg::BmEops> vec;

  while (fin.getline(line, 1024)) {
    int error = 0;
    double td[10];
    int sz = std::strlen(line);
    const char *str = line;
    for (int i = 0; i < 9; i++) {
      auto res = std::from_chars(skipws(str), line + sz, td[i]);
      if (res.ec != std::errc{})
        ++error;
      str = res.ptr;
    }
    if (!error) {
      int imjd = (int)td[0];
      double fsec = (td[0] - imjd) * 86400e0;
      const dso::MjdEpoch t(imjd, dso::FractionalSeconds{fsec});
      vec.emplace_back(t, &td[1]);
    }
  }

  if (!fin.eof()) {
    fprintf(stderr, "ERROR Failed parsing input file %s\n", fn);
    throw std::runtime_error("Failed parsing input file " + std::string(fn) +
                             "\n");
  }

  return vec;
}

const char *costg::basename(const char *fn) {
  const char *last_split = fn;
  const char *c = fn;
  while (*++c)
    if (*c == '/')
      last_split = c;
  return last_split + 1;
}
