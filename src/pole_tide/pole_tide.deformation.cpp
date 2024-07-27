#include "pole_tide.hpp"

dso::CartesianCrd
dso::PoleTide::deformation(const dso::MjdEpoch &t, double xp, double yp,
                           const dso::SphericalCrdConstView &rsta) noexcept {
  const auto m12 = pole_tide_details::mcoeffs(t, xp, yp);
  const double m1 = m12.m1; /* m1 in [arcsec] */
  const double m2 = m12.m2; /* m2 in [arcsec] */
  /* polar distance (co-latitude) */
  const double theta = dso::DPI / 2e0 - rsta.lat();
  /* trigonometric numbers */
  const double sl = std::sin(rsta.lon());
  const double cl = std::cos(rsta.lon());
  const double st = std::sin(theta);
  const double ct = std::cos(theta);
  const double s2t = 2e0 * st * ct;
  const double c2t = 2e0 * ct * ct - 1e0;
  /* IERS 2010, Eq. (26) */
  Eigen::Matrix<double, 3, 1> s;
  s(0) = -9e0 * c2t * (m1 * cl + m2 * sl);  /* Sθ along colatitude */
  s(1) = 9e0 * ct * (m1 * sl - m2 * cl);    /* Sλ along longitude */
  s(2) = -33e0 * s2t * (m1 * cl + m2 * sl); /* Sr radial */
  /* Rotation matrix to Cartesian coordinates; IERS 2010, Eq. (28) */
  Eigen::Matrix<double, 3, 3> R;
  R(0, 0) = ct * cl;
  R(1, 0) = ct * sl;
  R(2, 0) = -st;
  R(0, 1) = -sl;
  R(1, 1) = cl;
  R(2, 1) = 0e0;
  R(0, 2) = st * cl;
  R(1, 2) = st * sl;
  R(2, 2) = ct;
  /* apply rotation to transform to Cartesian */
  return dso::CartesianCrd(R * s);
}
