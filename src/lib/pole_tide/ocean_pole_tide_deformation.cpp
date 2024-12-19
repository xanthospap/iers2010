#include "pole_tide.hpp"
#include "geodesy/units.hpp"
#include "geodesy/core/crd_transformations.hpp"

dso::CartesianCrd
dso::OceanPoleTide::deformation(const MjdEpoch &t, double xp, double yp,
                                const dso::SphericalCrdConstView rsta,
                                const dso::OceanPoleTideDesaiCoeffs &coef,
                                double Re, double GM, double OmegaEarth,
                                double G, double ge) noexcept {

  /* density of sea water in [kgm^−3] */
  constexpr const double rhow = 1025e0;
  constexpr const double g2_real = 0.6870e0;
  constexpr const double g2_imag = 0.0036e0;

  const auto m12 = dso::pole_tide_details::mcoeffs(t, xp, yp);
  const double m1 = dso::sec2rad(m12.m1); /* m1 in [rad] */
  const double m2 = dso::sec2rad(m12.m2); /* m2 in [rad] */

  const double Hp = std::sqrt(4e0 * D2PI / 15) * (OmegaEarth * OmegaEarth) *
                    std::pow(Re, 4) / GM;
  const double K = 2e0 * D2PI * G * Re * rhow * Hp / 3e0 / ge;

    /* IERS 2010, Ch. 7.1.5, Eq. (29) */
  Eigen::Matrix<double, 3, 1> s;
  /* east */
  s(0) = (m1 * g2_real + m2 * g2_imag) * coef.eR +
         (m2 * g2_real - m1 * g2_imag) * coef.eI;
  /* north */
  s(1) = (m1 * g2_real + m2 * g2_imag) * coef.nR +
         (m2 * g2_real - m1 * g2_imag) * coef.nI;
  /* radial */
  s(2) = (m1 * g2_real + m2 * g2_imag) * coef.rR +
         (m2 * g2_real - m1 * g2_imag) * coef.rI;

  /* Rotation matrix to Cartesian coordinates (from spherical lon, lat) */
  const auto R = dso::geodetic2lvlh(rsta.lat(), rsta.lon());

  /* apply rotation to transform to Cartesian */
  return dso::CartesianCrd(K*(R * s));
}
