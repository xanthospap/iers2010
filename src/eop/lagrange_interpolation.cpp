#include "eop.hpp"
#include <datetime/dtdatetime.hpp>
#include <cstdio>

dso::EopSeries::EopInterpoationResult
dso::EopSeries::interpolate(const dso::MjdEpoch &t, dso::EopRecord &eop,
                            int order) const noexcept {
  /* find a matching interval */
  auto it = this->upper_bound(t);

  /* given date is prior to first date in EOP series */
  if (it == mvec.begin()) {
    /* copy the first entry in series and flag the return status */
    eop = *it;
    eop.t() = t;
    //fprintf(stderr, "--> prior to first eop record!\n");
    return EopInterpoationResult::OutOfBoundsPrior;
  }

  /* given date is later than last date in EOP series */
  if (it == mvec.end()) {
    --it;
    eop = *it;
    eop.t() = t;
    //fprintf(stderr, "--> after last eop record!\n");
    return EopInterpoationResult::OutOfBoundsLater;
  }

  /* check that we have enough data points (left & right) */
  //fprintf(stderr, "%.12f <= %.12f < %.12f\n", (it-1)->t().as_mjd(), t.as_mjd(), it->t().as_mjd());
  const int window = order + 1;
  int ndp = window / 2;
  //fprintf(stderr, "initial window = %d\n", window);
  if (it - ndp < mvec.begin() || it + (ndp - 1) >= mvec.end()) {
    /* not enough data points! try setting fewer data points, until 1 */
    --ndp;
    while ((ndp > 0) &&
           (it - ndp < mvec.begin() || it + (ndp - 1) >= mvec.end())) {
      --ndp;
    }
  }
  //fprintf(stderr, "reached window = %d\n", ndp*2);

  /* set all EOPs to 0 for output */
  eop.xp() = eop.yp() = eop.dut() = eop.lod() = eop.dX() = eop.dY() =
      eop.xp_rate() = eop.yp_rate() = 0e0;
  eop.t() = t;

  /* perform linear interpolation [it-1, it+1)*/
  if (ndp == 1) {
    const auto x0 = it-1;
    const auto x1 = it;
    const double dx1x0 =
        x1->t().diff<dso::DateTimeDifferenceType::FractionalDays>(x0->t());
    const double dx1x =
        x1->t().diff<dso::DateTimeDifferenceType::FractionalDays>(t);
    const double dxx0 =
        t.diff<dso::DateTimeDifferenceType::FractionalDays>(x0->t());
    const double dx0 = dx1x / dx1x0;
    const double dx1 = dxx0 / dx1x0;
    eop.xp()  += x0->xp()  * dx0 + x1-> xp() *dx1;
    eop.yp()  += x0->yp()  * dx0 + x1-> yp() *dx1;
    eop.dut() += x0->dut() * dx0 + x1-> dut()*dx1;
    eop.lod() += x0->lod() * dx0 + x1-> lod()*dx1;
    eop.dX()  += x0->dX()  * dx0 + x1-> dX() *dx1;
    eop.dY()  += x0->dY()  * dx0 + x1-> dY() *dx1;
    return EopInterpoationResult::Linear;
  }

  /* Lagrange interpolation with a total of ndp data points, from 
   * [it-ndp, it+ndp).
   */
  auto i = it - ndp;
  const auto stop = it + ndp;
  //fprintf(stderr, "start interpolation at %d and end at %d\n", (int)std::distance(mvec.begin(), i), (int)std::distance(mvec.begin(), stop));
  //double nom[ndp * 2]; /* Π_{k=0,k!=j}^{n} (x-x_k) */
  //double dom[ndp * 2]; /* Π_{k=0,k!=j}^{n} (x_j-x_k) */
  //double * __restrict__ nom = scr1();
  //double * __restrict__ dom = scr2();
  auto nom = work.begin();
  auto dom = work.begin() + MAX_POLY_INTERPOLATION_DEGREE+1;
  for (; i != stop; ++i) {
    /* j index runs through [0, 2ndp) */
    const int j = std::distance(it - ndp, i);
    /* x_j i.e. t_j */
    const auto xjt = i->t();
    nom[j] = 1e0;
    dom[j] = 1e0;
    auto k = it - ndp;
    for (; k != stop; ++k) {
      if (j != std::distance(it - ndp, k)) {
        const double n =
            t.diff<dso::DateTimeDifferenceType::FractionalDays>(k->t());
        const double d =
            xjt.diff<dso::DateTimeDifferenceType::FractionalDays>(k->t());
        nom[j] *= n;
        dom[j] *= d;
      }
    }
  }

  /* interpolate for different EOPs */
  i = it - ndp;
  for (; i != stop; ++i) {
    int j = std::distance(it - ndp, i);
    eop.xp() += i->xp() * nom[j] / dom[j];
    eop.yp() += i->yp() * nom[j] / dom[j];
    eop.dut() += (i->dut()-i->dat()) * nom[j] / dom[j];
    eop.lod() += i->lod() * nom[j] / dom[j];
    eop.dX() += i->dX() * nom[j] / dom[j];
    eop.dY() += i->dY() * nom[j] / dom[j];
  }

  eop.dut() += dso::dat(dso::modified_julian_day(t.tt2utc().imjd()));

  return (ndp * 2 == order - 1) ? EopInterpoationResult::PolyDegreeRequested
                                : EopInterpoationResult::PolyDegreeDescreased;
}
