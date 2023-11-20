#include "icgemio.hpp"
#include <charconv>

int dso::resolve_icgem_data_date_v1(const char *line, dso::Icgem::Datetime &t,
                                    const char *&ptr) noexcept {
  int ints[3];
  const char *str = line;
  int error = 0;

  /* Format: yyyymmdd */
  auto res = std::from_chars(str, str + 4, ints[0]);
  if (res.ec != std::errc{})
    ++error;
  str = res.ptr;

  res = std::from_chars(str, str + 2, ints[1]);
  if (res.ec != std::errc{})
    ++error;
  str = res.ptr;

  res = std::from_chars(str, str + 2, ints[2]);
  if (res.ec != std::errc{})
    ++error;

  /* check for error(s) */
  if (error) {
    fprintf(stderr,
            "[ERROR] Failed resolving ICGEM v2 date: %.13s (traceback: %s)\n",
            line, __func__);
    return error;
  }

  /* construct datetime */
  t = dso::Icgem::Datetime(dso::year(ints[0]), dso::month(ints[1]),
                           dso::day_of_month(ints[2]), dso::hours(0),
                           dso::minutes(0), dso::nanoseconds(0));

  /* set pointer to first, unresolved char */
  ptr = res.ptr;

  /* all done */
  return 0;
}
