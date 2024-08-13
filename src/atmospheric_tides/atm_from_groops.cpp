#include <vector>
#include <string>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <cstring>
#include "doodson.hpp"
#include <charconv>

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

struct GroopsTideModelFileName {
  char _model[24]={"\0"};
  dso::DoodsonConstituent _doodson;
  char _constituentName[8]={"\0"};
  char _sinCos;

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

int dso::atm_tide_from_groops(const char *file_list, const char *doodson,
                              const char *admittance,
                              const char *data_dir) noexcept {
}
