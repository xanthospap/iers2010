#include "doodson.hpp"
#include <charconv>
#include <cstdio>
#include <stdexcept>

dso::DoodsonConstituent dso::resolve_iers10_doodson_string(const char *cstr) {
  int ints[6];
  int error = 0;

  const char *str = cstr;
  while (*str && *str == ' ')
    ++str;

  ints[0] = *str++ - '0';
  ints[1] = *str++ - '0';
  ints[2] = *str++ - '0';
  if (*str != '.')
    error = 1;
  
  ++str;
  ints[3] = *str++ - '0';
  ints[4] = *str++ - '0';
  ints[5] = *str++ - '0';

  for (int i = 0; i < 6; i++)
    if (ints[i] < 0 || ints[i] > 9)
      error += 1;

  if (error) {
    fprintf(stderr,
            "[ERROR] Failed resolving Doodson number from string: %s "
            "(traceback: %s)\n",
            cstr, __func__);
    throw std::runtime_error("[ERROR] Failed resolving Doodson number\n");
  }

  return dso::DoodsonConstituent(ints);
}
