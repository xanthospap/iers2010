#include "earth_rotation.hpp"
#include "iau.hpp"

Eigen::Vector3d dso::c2i06a(const dso::MjdEpoch& tt, const dso::EopRecord& eops, Eigen::Quaterniond& qc2tirs, Eigen::Quaterniond& qtirs2i) noexcept
{
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

  /* quaternion GCRS-to-TIRS */
  qc2tirs = dso::detail::c2tirs(dso::era00(tt.tt2ut1(eops.dut())), s, d, e);

  /* quaternion TIRS-to-ITRS */
  qtirs2i = dso::detail::tirs2i(dso::sec2rad(eops.xp()), dso::sec2rad(eops.yp()), sp00(tt));

  /* rotation vector (ω_x, ω_y, ω_z) */
  return dso::earth_rotation_axis(eops.lod());
}
