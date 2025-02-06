#include "aod1b.hpp"
#include "datetime/datetime_write.hpp"
#include <cstdio>
#include <cstring>
#include <charconv>

const char *dso::Aod1bNonTidalProductNaming::filename(
    const dso::ymd_date &t, int rl, char *buf, char satid, bool is_compressed) noexcept {
  /* first 6 chars */
  std::strcpy(buf, "AOD1B_");
  /* adding YYYY-MM-DD, i.e. another 10 chars */
  dso::to_char<dso::YMDFormat::YYYYMMDD>(t, buf + 6, '-');
  /* adding satellite id, i.e. the part '_S' (another 2 chars) */
  buf[16] = '_';
  buf[17] = satid;
  /* adding product '_RL.EXT' (another 7 chars) */
  std::sprintf(buf + 18, "_%02d.asc", rl);
  /* add compression (if needed) */
  if (is_compressed) 
    std::strcpy(buf+27, ".gz");
  /* all done */
  return buf;
}

int dso::Aod1bNonTidalProductNaming::resolve_rl(const char *fn) noexcept {
  /* Find the last occurrence of '/' to get the filename */
  const char *basename = std::strrchr(fn, '/');
  ++basename;
  int rl;
  const auto j = std::from_chars(basename + 19, basename + 21, rl);
  if (j.ec != std::errc{}) {
    fprintf(
        stderr,
        "[ERROR] Failed resolving RL from AOD1B filename %s (traceback: %s)\n",
        basename, __func__);
    return -100;
  }
  return rl;
}
