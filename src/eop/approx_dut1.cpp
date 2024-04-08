#include "eop.hpp"

dso::EopSeries::EopInterpolationResult
dso::EopSeries::approx_dut1(const dso::MjdEpoch &t,
                            double &dut1) const noexcept {
  /* find a matching interval (try for last used iterator first) */
  auto it = last_it;
  {
    if (last_it != mvec.begin() && last_it != mvec.end()) {
      if (!((last_it - 1)->t() >= t && last_it->t() < t)) {
        it = this->upper_bound(t);
      }
    } else {
      it = this->upper_bound(t);
    }
  }
  last_it = it;

  /* given date is prior to first date in EOP series */
  if (it == mvec.begin()) {
    /* copy the first entry in series and flag the return status */
    dut1 = it->dut();
    return EopInterpolationResult::OutOfBoundsPrior;
  }

  /* given date is later than last date in EOP series */
  if (it == mvec.end()) {
    --it;
    dut1 = it->dut();
    return EopInterpolationResult::OutOfBoundsLater;
  }

  /* perform linear interpolation [it-1, it+1) */
  const auto x0 = it - 1;
  const auto x1 = it;
  const double dx1x0 =
      x1->t().diff<dso::DateTimeDifferenceType::FractionalDays>(x0->t());
  const double dx1x =
      x1->t().diff<dso::DateTimeDifferenceType::FractionalDays>(t);
  const double dxx0 =
      t.diff<dso::DateTimeDifferenceType::FractionalDays>(x0->t());
  const double dx0 = dx1x / dx1x0;
  const double dx1 = dxx0 / dx1x0;
  /* continuous (UT1-TAI)j <- (UT1-UTC)j - ΔATj */
  dut1 = (x0->dut() - x0->dat()) * dx0 + (x1->dut() - x1->dat()) * dx1;
  dut1 += dso::dat(dso::modified_julian_day(t.tt2utc().imjd()));

  return EopInterpolationResult::Linear;
}
