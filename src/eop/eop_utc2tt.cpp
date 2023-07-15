#include "eop.hpp"
#include <algorithm>

void dso::EopLookUpTable::utc2tt() noexcept {
  std::transform(tvec.cbegin(), tvec.cend(), tvec.begin(),
                 [](const dso::TwoPartDate &utc) { return utc.utc2tt(); });
  return;
}
