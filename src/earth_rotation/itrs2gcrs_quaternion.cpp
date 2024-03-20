#include "earth_rotation.hpp"
#include "iau.hpp"
#include "geodesy/units.hpp"
#include <stdexcept>
#include <cstring>

Eigen::Quaterniond
dso::EarthRotation::itrs2gcrs_quaternion(const dso::MjdEpoch &tt,
                                         double *fargs14) {
  /* interpolate EOPs for the given epoch (store in meop) */
  if (eops_at(tt)) {
    throw std::runtime_error(
        "[ERROR] Faled to interpolate EOP data for given epoch\n");
  }

  /* Eart Rotation Angle (need UT1, and ΔUT) */
  const double era = era00(tt.tt2ut1(meop.dut()));

  /* CIP components */
  double fargs[14];
  double xcip, ycip;
  dso::xycip06a(tt, xcip, ycip, fargs);
  /* Note: Xcip = X(IAU 2006/2000A) + δX */
  xcip += dso::sec2rad(meop.dX());
  ycip += dso::sec2rad(meop.dY());

  /* CIO and TIO locator */
  const double s = dso::s06(tt, xcip, ycip, fargs14);
  const double sp = dso::sp00(tt);

  /* shall we return the fundamental arguemnts ? */
  if (fargs14) std::memcpy(fargs14, fargs, sizeof(double)*14);

  /* compute rotation quaternion */
  return dso::detail::itrs2gcrs_quaternion(
      era, s, sp, xcip, ycip, dso::sec2rad(meop.xp()), dso::sec2rad(meop.yp()));
}
