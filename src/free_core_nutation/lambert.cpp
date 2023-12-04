#include "fcn.hpp"
#include "geodesy/geoconst.hpp"
#include <algorithm>

dso::fcn::FcnResult
dso::lambert_fcn(const dso::MjdEpoch &t,
                 const std::vector<fcn::LambertCoefficients> &lvec) noexcept {
  /* Mean prediction error in [microas / day] */
  constexpr const double mpe = 0.3e0;

  /* FCN parameter period in days */
  constexpr const double per = -430.21e0;

  /* phase in [rad] at given date; time is fractional days since J2000 */
  const double phi =
      (dso::D2PI / per) * (t.diff<dso::DateTimeDifferenceType::FractionalDays>(
                              dso::MjdEpoch::j2000_mjd()));

  /* find an entry at the coefficient vector, for which t0 <= t < t1 */
  auto rit = lvec.rbegin() + 1;
  for (; rit != lvec.rend(); ++rit) {
    if (t >= rit->t && t < (rit - 1)->t)
      break;
  }

  double axc[2], axs[2], ayc[2], ays[2];
  /* interval found */
  if (rit != lvec.rend()) {
    const auto prev = *rit;
    const auto next = *(rit - 1);
    const double dt =
        t.diff<dso::DateTimeDifferenceType::FractionalDays>(prev.t);
    const double dT =
        next.t.diff<dso::DateTimeDifferenceType::FractionalDays>(prev.t);
    /* CIP coordinates */
    const double daxc = next.xc() - prev.xc();
    const double daxs = next.xs() - prev.xs();
    const double dayc = next.yc() - prev.yc();
    const double days = next.ys() - prev.ys();
    axc[0] = prev.xc() + (daxc / dT) * dt;
    axs[0] = prev.xs() + (daxs / dT) * dt;
    ayc[0] = prev.yc() + (dayc / dT) * dt;
    ays[0] = prev.ys() + (days / dT) * dt;
    /* sigma values */
    const double sdaxc = next.sx() - prev.sx();
    const double sdayc = next.sy() - prev.sy();
    axc[1] = std::abs(prev.sx() + (sdaxc / dT) * dt);
    axs[1] = std::abs(prev.sx() + (sdaxc / dT) * dt);
    ayc[1] = std::abs(prev.sy() + (sdayc / dT) * dt);
    ays[1] = std::abs(prev.sy() + (sdayc / dT) * dt);
  } else {
    /* date is later than the most recent one */
    if (t >= lvec.rbegin()->t) {
      const auto it = *lvec.rbegin();
      const double dt =
          t.diff<dso::DateTimeDifferenceType::FractionalDays>(it.t);
      axc[0] = it.xc();
      axs[0] = it.xs();
      ayc[0] = it.yc();
      ays[0] = it.ys();
      axc[1] = it.sx() + mpe * dt;
      axs[1] = it.sx() + mpe * dt;
      ayc[1] = it.sy() + mpe * dt;
      ays[1] = it.sy() + mpe * dt;
    } else {
      /* date is earlier than the oldest date in vector */
#ifdef DEBUG
      assert(t < lvec.begin()->t);
#endif
      const auto it = *lvec.begin();
      const double dt =
          t.diff<dso::DateTimeDifferenceType::FractionalDays>(it.t);
      axc[0] = it.xc();
      axs[0] = it.xs();
      ayc[0] = it.yc();
      ays[0] = it.ys();
      axc[1] = it.sx() + mpe * dt;
      axs[1] = it.sx() + mpe * dt;
      ayc[1] = it.sy() + mpe * dt;
      ays[1] = it.sy() + mpe * dt;
    }
  }

  /* Computation of X and Y */
  const double cosp = std::cos(phi);
  const double sinp = std::sin(phi);
  const double X = axc[0] * cosp - axs[0] * sinp;
  const double Y = ayc[0] * cosp - ays[0] * sinp;

  /* Computation of the uncertainties */
  const double dX = axc[1] + axs[1];
  const double dY = ayc[1] + ays[1];

  /* accumulate and return result */
  return dso::fcn::FcnResult(X, Y, dX, dY);
}
