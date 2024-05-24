#include "doodson.hpp"
#include <cstring>

char *
dso::DoodsonConstituent::str(char *buf,
                             bool use_5s_convention = false) const noexcept {
  const int add = use_5s_convention ? 5 : 0;
  std::sprintf(buf, "%d%d%d.%d%d%d", iar[0], iar[1] + add, iar[2] + add,
               iar[3] + add, iar[4] + add, iar[5] + add);
  return buf;
}
