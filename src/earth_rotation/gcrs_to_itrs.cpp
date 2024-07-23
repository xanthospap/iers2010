#include "earth_rotation.hpp"
#include "iau.hpp"

Eigen::Quaterniond dso::itrs2gcrs_quaternion(const dso::MjdEpoch &tt,
                                             const dso::EopRecord &eops,
                                             double *fargs) noexcept {
  /* CIP and CIO, IAU 2006/2000A (units: [rad]). */
  double xcip, ycip;
  dso::xycip06a(tt, xcip, ycip, fargs);
  const double s = dso::s06(tt, xcip, ycip);

  /* Earth rotation Angle [rad] */
  const double era = dso::era00(tt.tt2ut1(eops.dut()));

  /* TIO locator s' [rad] */
  const double sp = dso::sp00(tt);

  return dso::detail::gcrs2itrs_quaternion(era, s, sp, xcip, ycip, eops.xp(),
                                           eops.yp());
}
