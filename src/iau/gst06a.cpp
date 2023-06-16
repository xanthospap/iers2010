#include "iau.hpp"

double iers2010::sofa::gst06a(const dso::TwoPartDate &mjd_ut1,
                     const dso::TwoPartDate &mjd_tt) noexcept {
  /* Classical nutation x precession x bias matrix, IAU 2000A. */
  const auto rnpb = pnm06a(mjd_tt);
  /* Greenwich apparent sidereal time. */
  return gst06(mjd_ut1, mjd_tt, rnpb);
}
