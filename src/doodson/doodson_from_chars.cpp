#include "doodson.hpp"
#include <cstdio>
#include <stdexcept>

dso::DoodsonConstituent dso::DoodsonConstituent::from_chars(const char *cstr) {
  int ints[6] = {
    dso::DoodsonConstituent::char2int(cstr[0]),
    dso::DoodsonConstituent::char2int(cstr[1])-5,
    dso::DoodsonConstituent::char2int(cstr[2])-5,
    dso::DoodsonConstituent::char2int(cstr[4])-5,
    dso::DoodsonConstituent::char2int(cstr[5])-5,
    dso::DoodsonConstituent::char2int(cstr[6])-5};

  if (cstr[3]!='.') {
    fprintf(stderr,
            "[ERROR] Failed resolving Doodson number from string: %s "
            "(traceback: %s)\n",
            cstr, __func__);
    throw std::runtime_error("[ERROR] Failed resolving Doodson number\n");
  }

  return dso::DoodsonConstituent(ints);
}
