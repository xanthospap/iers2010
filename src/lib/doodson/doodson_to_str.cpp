#include "doodson.hpp"
#include <cstring>

char *dso::DoodsonConstituent::str(char *buf,
                                   bool use_5s_convention) const noexcept {
  const int add = use_5s_convention ? 5 : 0;
  std::sprintf(buf, "%c%c%c.%c%c%c", dso::DoodsonConstituent::int2char(iar[0]),
               dso::DoodsonConstituent::int2char(iar[1] + add),
               dso::DoodsonConstituent::int2char(iar[2] + add),
               dso::DoodsonConstituent::int2char(iar[3] + add),
               dso::DoodsonConstituent::int2char(iar[4] + add),
               dso::DoodsonConstituent::int2char(iar[5] + add));
  return buf;
}
