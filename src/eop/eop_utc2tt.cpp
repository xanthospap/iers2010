#include "eop.hpp"
#include <algorithm>
#include <datetime/dtcalendar.hpp>

void dso::EopLookUpTable::utc2tt() noexcept {
  std::transform(t.cbegin(), t.cend(), t.begin(),
                 [](const dso::TwoPartDate &utc) { return utc.utc2tt(); });
  return;
}
