#include "iau.hpp"

double iers2010::sofa::gst06(const dso::TwoPartDate &mjd_ut1,
                             const dso::TwoPartDate &mjd_tt,
                             const Eigen::Matrix<double, 3, 3> &rnpb) noexcept {
  // Extract CIP coordinates.
  // iauBpn2xy(rnpb, &x, &y);
  const double x = rnpb(2, 0); //[2][0]
  const double y = rnpb(2, 1); //[2][1]

  // The CIO locator, s.
  const double s = iers2010::sofa::s06(mjd_tt, x, y);

  // Greenwich apparent sidereal time.
  const double era = iers2010::sofa::era00(mjd_ut1);
  const double eors = iers2010::sofa::eors(rnpb, s);
  return dso::anp(era - eors);
}
