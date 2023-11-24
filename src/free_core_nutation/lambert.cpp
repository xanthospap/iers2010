#include "fcn.hpp"
#include "geodesy/geoconst.hpp"
#include <algorithm>

int dso::lambert_fcn(const dso::TwoPartDate &t, const std::vector<fcn::LambertCoefficients> &lvec) noexcept {
  /* Mean prediction error in [microas / day] */
  constexpr const double mpe = 0.3e0;

  /* FCN parameter period in days */
  constexpr const double per = -430.21e0;

  /* phase in [rad] at given date; time is fractional days since J2000 */
  const double phi =
      (dso::D2PI / per) * (t.diff<dso::DateTimeDifferenceType::FractionalDays>(
                              dso::TwoPartDate::j2000_mjd()));

  /* find an entry at the coefficient vector, for which t0 <= t < t1 */
  const auto rit = lvec.rbegin()+1;
  for (; rit != lvec.rend(); ++rit) {
    if (t >= rit->t && t < (rit-1)->t)
      break;
  }

  double daxc,daxs,dayc,days,axc,axs,ayc,ays;
  /* interval found */
  if (rit != lvec.rend()) {
    const auto prev = *rit;
    const auto next = *(rit-1);
    const double dt = t.diff<dso::DateTimeDifferenceType::FractionalDays>(prev.t);
    const double dT = next.t.diff<dso::DateTimeDifferenceType::FractionalDays>(prev.t);
    daxc = next.xc() - prev.xc();
    daxs = next.xs() - prev.xs();
    dayc = next.yc() - prev.yc();
    days = next.ys() - prev.ys();
    axc = prev.xc() + (daxc / dT) * dt;
    axs = prev.xs() + (daxs / dT) * dt;
    ayc = prev.yc() + (dayc / dT) * dt;
    ays = prev.ys() + (days / dT) * dt;
  } else {
    /* date is later than the most recent one */
    if (t>= lvec.rbegin()->t) {
         const auto it = *lvec.rbegin();
         const double dt =  t.diff<dso::DateTimeDifferenceType::FractionalDays>(it.t); 
         axc=it.sx()+mpe*dt;
         axs=it.sx()+mpe*dt;
         ayc=it.sy()+mpe*dt;
         ays=it.sy()+mpe*dt;
    } else {
      /* date is earlier than the oldest date in vector */
#ifdef DEBUG
      assert(t< lvec.begin()->t);
#endif
         const auto it = *lvec.begin();
         const double dt =  t.diff<dso::DateTimeDifferenceType::FractionalDays>(it.t); 
         axc=it.sx()+mpe*dt;
         axs=it.sx()+mpe*dt;
         ayc=it.sy()+mpe*dt;
         ays=it.sy()+mpe*dt;
     }
  }

}
