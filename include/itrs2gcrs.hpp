/** @file
 * Define a class to handle the ITRS to GCRS (and the inverse) transformation,
 * both for position and velocity vector
 */
#ifdef TODO
#ifndef __DSO_ITRS_TO_GCRS_ROTATION_HPP__
#define __DSO_ITRS_TO_GCRS_ROTATION_HPP__

#include "eop.hpp"
#include "iau.hpp"

namespace dso {

class RotItrsGcrs {
private:
  /** Last date requested in TT */
  MjdEpoch mtt;
  /** Last EOPs interpolated, valid at mtt */
  EopRecord meops;
  /** An EopSeries instance, i.e. EOP data to compute rotation */
  EopSeries meopsv;
  /** fundamental arguments, lunisolar + planetary at mtt */
  double fargs[14];

  /** @brief Use the EopSeries to interpolate EOP values at the given TT */
  void get_eops(const MjdEpoch &tt) noexcept {
    meopsv.interpolate(tt, meops);
    mtt = tt;
  }

public:
  /** @brief Get the ITRS to GCRS quaternion */
  Eigen::Quaterniond get_quaternion(const MjdEpoch &tt) noexcept {
    if (tt != mtt) {
      meopsv.interpolate(tt, meops);
      /* fundamental arguments, lunisolar and planetary */
      fargs[0] = fal03(mtt);
      fargs[1] = falp03(mtt);
      fargs[2] = faf03(mtt);
      fargs[3] = fad03(mtt);
      fargs[4] = faom03(mtt);
      fargs[5] = fame03(mtt);
      fargs[6] = fave03(mtt);
      fargs[7] = fae03(mtt);
      fargs[8] = fama03(mtt);
      fargs[9] = faju03(mtt);
      fargs[10] = fasa03(mtt);
      fargs[11] = faur03(mtt);
      fargs[12] = fane03(mtt);
      fargs[13] = fapa03(mtt);
      /* get (X,Y) CIP and s */
      double xcip, ycip;
      xycip06a(tt, xcip, ycip, fargs);
      const double s = s06(tt, xcip, ycip);
      /* add (interpolated) observed correction to (X,Y) CIP */
      xcip += dso::sec2rad(meops.dX());
      ycip += dso::sec2rad(meops.dY());
      /* get the ERA */
      const double era = era00(tt.tt2ut1(meops.dut()));
    }
    return itrs2gcrs_quaternion(era, s, sp00(tt), xcip, ycip, meops.xp(),
                                meops.yp());
  }

}; /* RotItrsGcrs */

} /* namespace dso */

#endif
#endif
