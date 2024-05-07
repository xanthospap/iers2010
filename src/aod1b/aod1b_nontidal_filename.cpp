#include "aod1b.hpp"
#include "datetime/datetime_read.hpp"
#include "datetime/datetime_write.hpp"
#include <cstring>
#include <stdexcept>
#include <filesystem>

namespace {
  inline
  const char *mbasename(const char *path) noexcept {
    const auto sep = std::filesystem::path::preferred_separator;
    int sz = std::strlen(path);
    const char *c = path + sz - 1;
    while ((c != path) && (*c != sep)) --c;
    return (c==path)?(path):(c+1);
  }
} /* unnamed namespace */

dso::Aod1bIn::Aod1bNonTidalFilenameStruct::Aod1bNonTidalFilenameStruct(
    const char *fn_aod1b) {
  /* first off, we need to strip the path */
  const char *fn = mbasename(fn_aod1b);
  if ((std::strlen(fn) < 25) || std::strncmp(fn, "AOD1B_", 6)) {
    fprintf(stderr,
            "[ERROR] Failed resolving AOD1B filename %s (traceback: %s)\n", fn,
            __func__);
    throw std::runtime_error("[ERROR] Failed resolving AOD1B filename\n");
  }
  /* read date */
  const char *end;
  const auto ymd =
      dso::ReadInDate<dso::YMDFormat::YYYYMMDD>::read(fn + 6, &end);
  mt = dso::Datetime<dso::seconds>(ymd);
  /* read satellite (should be 'X') */
  msatid = *(++end);
  /* read revision/version */
  ++end;
  mrevision[0] = *(++end);
  mrevision[1] = *(++end);
  /* read extension */
  assert(*(++end) == '.');
  std::memcpy(mext, end + 1, 3);
}

const char *
dso::Aod1bIn::Aod1bNonTidalFilenameStruct::make(char *buf) const noexcept {
  std::memcpy(buf, "AOD1B_", 6);
  dso::to_char<dso::YMDFormat::YYYYMMDD>(mt.as_ymd(), buf + 6, '-');
  std::sprintf(buf + 16, "_%c_%.2s.%.3s", msatid, mrevision, mext);
  return buf;
}
