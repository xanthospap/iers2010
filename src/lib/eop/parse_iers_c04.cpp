#include "eop.hpp"

int dso::parse_iers_C04(const char *c04fn, const MjdEpoch &start_tt,
                        const dso::MjdEpoch &end_tt,
                        dso::EopSeries &eops) noexcept {
  dso::details::IersEopFormat type;

  if (dso::details::choose_c04_series(c04fn, type))
    return 1;

  int error;
  if (type == dso::details::IersEopFormat::C0420)
    error = dso::details::parse_iers_C0420(c04fn, start_tt, end_tt, eops);
  else
    error = dso::details::parse_iers_C0414(c04fn, start_tt, end_tt, eops);

  eops.reset_vec_iterator();
  return error;
}
