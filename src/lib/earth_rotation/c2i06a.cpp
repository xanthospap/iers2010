#include "earth_rotation.hpp"
#include "iau.hpp"

Eigen::Quaterniond dso::c2i06a(const dso::MjdEpoch &tt,
                               const dso::EopRecord &eops) noexcept {
  using namespace Eigen;

  double fargs[14];
  double Xcip, Ycip;

  /* compute (X,Y) CIP and fundamental arguments */
  dso::xycip06a(tt, Xcip, Ycip, fargs);

  /* use fundamental arguments to compute s */
  const double s = dso::s06(tt, Xcip, Ycip, fargs);

  /* apply CIP corrections */
  Xcip += dso::sec2rad(eops.dX());
  Ycip += dso::sec2rad(eops.dY());

  /* spherical crd for CIP */
  double d, e;
  dso::detail::xycip2spherical(Xcip, Ycip, d, e);

  /* accumulated rotation */
  return dso::detail::c2i(dso::era00(tt.tt2ut1(eops.dut())), s, dso::sp00(tt),
                          d, e, dso::sec2rad(eops.xp()),
                          dso::sec2rad(eops.yp()));
}