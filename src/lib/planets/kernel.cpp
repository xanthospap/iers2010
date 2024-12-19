#include "planets.hpp"
#include <cstring>

namespace {
constexpr const int FILLEN = 128;
constexpr const int TYPLEN = 32;
constexpr const int SRCLEN = 128;
}

int dso::cspice::load_if_unloaded(const char *kernel) noexcept {
  SpiceInt count, which, handle;
  SpiceChar file[FILLEN];
  SpiceChar filtyp[TYPLEN];
  SpiceChar source[SRCLEN];
  SpiceBoolean found;

  /* get number of all ascii kernels loaded ... (to find if we already have
   * the LSK kernel loaded)
   */
  ktotal_c("ALL", &count);

  if (count) {
    for (which = 0; which < count; which++) {
      kdata_c(which, "ALL", FILLEN, TYPLEN, SRCLEN, file, filtyp, source,
              &handle, &found);
      if (!std::strncmp(file, kernel, std::strlen(kernel)))
        return 0;
    }
  }

  furnsh_c(kernel);
  return 0;
}
