#include <vector>
#include <string>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <charconv>
#include "eigen3/Eigen/Eigen"
#include "doodson.hpp"
#include "icgemio.hpp"
#include "atmospheric_tides.hpp"

namespace {
/* given a file that contains two-column string, seperated by whitespace, 
 * resolve them and return them in a vector of strings. File example:
 * TiME22_ATM_162.556_pi1_cos.gfc TiME22_ATM_162.556_pi1_sin.gfc 
 * TiME22_ATM_163.555_p1_cos.gfc TiME22_ATM_163.555_p1_sin.gfc 
 * TiME22_ATM_164.556_s1_cos.gfc TiME22_ATM_164.556_s1_sin.gfc 
 * TiME22_ATM_165.555_k1_cos.gfc TiME22_ATM_165.555_k1_sin.gfc 
 * TiME22_ATM_166.554_psi1_cos.gfc TiME22_ATM_166.554_psi1_sin.gfc 
 * TiME22_ATM_255.555_m2_cos.gfc TiME22_ATM_255.555_m2_sin.gfc 
 * TiME22_ATM_272.556_t2_cos.gfc TiME22_ATM_272.556_t2_sin.gfc 
 * 
 * returns a vector<string> with elements:
 * TiME22_ATM_162.556_pi1_cos.gfc
 * TiME22_ATM_162.556_pi1_sin.gfc
 * TiME22_ATM_163.555_p1_cos.gfc
 * TiME22_ATM_163.555_p1_sin.gfc
 * TiME22_ATM_164.556_s1_cos.gfc
 * TiME22_ATM_164.556_s1_sin.gfc
 * TiME22_ATM_165.555_k1_cos.gfc
 * TiME22_ATM_165.555_k1_sin.gfc
 * TiME22_ATM_166.554_psi1_cos.gfc
 * TiME22_ATM_166.554_psi1_sin.gfc
 * TiME22_ATM_255.555_m2_cos.gfc
 * [...]
 */
std::vector<std::string> resolve_file_list(const char *fn) noexcept {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening file %s (traceback: %s)\n", fn,
            __func__);
    return std::vector<std::string>{};
  }

  std::vector<std::string> fnvec;
  constexpr const int LSZ = 512;
  char line[LSZ];
  while (fin.getline(line, LSZ)) {
    const char *str1 = line;
    while (*str1 && *str1==' ') ++str1;
    const char *end1 = str1;
    while (*end1 && *end1!=' ') ++end1;
    const char *str2 = end1;
    while (*str2 && *str2==' ') ++str2;
    const char *end2 = str2;
    while (*end2 && *end2!=' ') ++end2;
    std::string sline(line);
    fnvec.emplace_back(sline.substr(str1-line, end1-str1));
    fnvec.emplace_back(sline.substr(str2-line, end2-str2));
  }
  return fnvec;
}

const char *skipws(const char *line) noexcept {
  while (*line && *line == ' ') ++line;
  return line;
}

std::vector<std::array<int,6>> resolve_doodson(const char *fn) noexcept {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening file %s (traceback: %s)\n", fn,
            __func__);
    return std::vector<std::array<int,6>>{};
  }

  std::vector<std::array<int,6>> dvec;
  constexpr const int LSZ = 124;
  char line[LSZ];
  int error = 0;
  while (fin.getline(line, LSZ) && (!error)) {
    std::array<int,6> ar;
    auto sz = std::strlen(line);
    const char *str = line;
    for (int i=0; i<6; i++) {
      auto res = std::from_chars(skipws(str), str+sz, ar[i]);
      error += (res.ec != std::errc{});
      str = res.ptr;
    }
    dvec.emplace_back(ar);
  }

  if (error) {
    fprintf(stderr, "[ERROR] Failed resolving Doodson numbers from file %s (traceback: %s)\n", fn, __func__);
    return std::vector<std::array<int,6>>{};
  }

  return dvec;
}

int parse_admittance(const char *fn, Eigen::MatrixXd &mat) noexcept {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening file %s (traceback: %s)\n", fn,
            __func__);
    return 1;
  }

  const int rows = mat.rows();
  const int cols = mat.cols();

  constexpr const int LSZ = 1024;
  char line[LSZ];
  int error = 0;
  int row = 0;

  while (fin.getline(line, LSZ) && (!error)) {
    if (row > rows) {
      fprintf(stderr,
              "[ERROR] More rows than expected in admittance file %s "
              "(traceback: %s)\n",
              fn, __func__);
      return 1;
    }
    auto sz = std::strlen(line);
    const char *str = line;
    for (int col = 0; col < cols; col++) {
      auto res = std::from_chars(skipws(str), str + sz, mat(row, col));
      error += (res.ec != std::errc{});
      str = res.ptr;
    }
    ++row;
  }

  if (error) {
    fprintf(
        stderr,
        "[ERROR] Failed resolving admittance from file %s (traceback: %s)\n",
        fn, __func__);
    return 1;
  }

  return 0;
}

struct GroopsTideModelFileName {
  char _model[24]={"\0"};
  dso::DoodsonConstituent _doodson;
  char _constituentName[8]={"\0"};
  char _sinCos; /* 'c' or 's' */

  GroopsTideModelFileName(const char *fn) {
    int error = 0;
    /* pattern is <model name>_<doodson>_<name>_<cos/sin>.gfc */
    const int sz = std::strlen(fn);
    /* start backwards; check extension fn[sz-4:sz] == .gfc */
    if (std::strncmp(fn+sz-4, ".gfc", 4)) ++error;
    
    /* check sin/cos fn[sz-7] == sin/cos */
    {
      if (std::strncmp(fn + sz - 7, "sin", 3)) {
        if (std::strncmp(fn + sz - 7, "cos", 3)) {
          ++error;
        } else {
          _sinCos = 'c';
        }
      } else {
        _sinCos = 's';
      }
    }

    /* get name of harmonic */
    int lastchar = sz - 9; /* should be just be before '_[sincos]' */
    {
      const char *str = fn + lastchar;
      while (*str && *str != '_')
        --str;
      std::memcpy(_constituentName, str+1, fn + sz - 9 - str);
      // printf("\tconstituent: [%s]\n", _constituentName);
      lastchar = str - fn;
    }

    /* get Doodson */
    {
      const char *str = fn + lastchar - 1;
      while (*str && *str != '_')
        --str;
      //std::memcpy(buf, str+1, fn+lastchar-1-str);
      _doodson = dso::DoodsonConstituent::from_chars(str+1);
      lastchar = str-fn;
    }

    /* the rest is the model name */
    std::memcpy(_model, fn, lastchar);
  }
}; /* GroopsTideModelFileName */
} /* unnamed namespace */

int atm_tide_from_groops(const char *file_list, const char *doodson,
                              const char *admittance,
                              const char *data_dir, 
                              int max_degree, int max_order) noexcept {

  /* parse file list from <model>_001_fileList.txt */
  const auto flvec = resolve_file_list(file_list);
  if (flvec.size() < 1) {
    fprintf(stderr,
            "[ERROR] Failed parsing file list from file %s (traceback: %s)\n",
            file_list, __func__);
    return 1;
  }

  if (flvec.size() / 2) {
    fprintf(stderr,
            "[ERROR] Expected even number of files in input file %s, found %ld "
            "(traceback: %s)\n",
            file_list, flvec.size(), __func__);
    return 1;
  }

  /* tidal lines in tidal atlas (i.e. k) */
  int k = flvec.size() / 2;

  /* parse Doodson numbers; number of tidal waves (doodson numbers) is f */
  const auto ddvec = resolve_doodson(doodson);
  int f = ddvec.size();

  /* admittance matrix A of size: (k,f) */
  Eigen::MatrixXd A(k, f);
  if (parse_admittance(admittance, A)) {
    fprintf(stderr,
            "[ERROR] Failed parsing admittance matrix from file %s (traceback: "
            "%s)\n",
            admittance, __func__);
    return 1;
  }

  /* Create a new (empty) AtmosphericTide instance */
  dso::AtmosphericTide atm;

  /* iterate through input files (waves) and read/parse coefficients for each 
   * of the k waves of the atlas
   */
  for (const auto &wave_fn : flvec) {
    /* resolve filename */
    GroopsTideModelFileName wave_info(wave_fn.c_str());
    /* check if we already have the tidal wave (in AtmosphericTide instance) */
    auto it = std::find_if(atm.wave_vector().begin(), atm.wave_vector().end(),
                        [=](const dso::detail::AtmosphericTidalWave &w) {
                          return w.wave().doodson() == wave_info._doodson;
                        });
    if (it == atm.wave_vector().end()) {
      /* new tidal wave */
      it = atm.append_wave(
          dso::detail::TidalConstituentsArrayEntry(wave_info._doodson, 0e0, 0e0,
                                                   wave_info._constituentName),
          max_degree, max_order);
    } else {
      /* tidal wave already exists */
      ;
    }
    /* an Icgem instance to read data from (Stokes coefficients) */
    dso::Icgem icgem(wave_fn.c_str());
    /* read coefficients to sin/cos part */
    dso::StokesCoeffs *cs = &(it->stokes_cos());
        //(wave_info._sinCos == 's') ? &(it->stokes_sin()) : &(it->stokes_cos());
    if (icgem.parse_data(max_degree, max_order, dso::Icgem::Datetime::min(), *cs)) {
      fprintf(stderr,
              "[ERROR] Failed parsing Stokes coefficients of type \'%c\' from "
              "file %s (traceback: %s)\n",
              wave_info._sinCos, wave_fn.c_str(), __func__);
      return 1;
    }
  }
}
