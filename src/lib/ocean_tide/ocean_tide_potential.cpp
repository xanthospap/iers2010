#include "iau.hpp"
#include "ocean_tide.hpp"

int dso::OceanTide::stokes_coeffs(const dso::MjdEpoch& mjdtt,
    const dso::MjdEpoch& mjdut1,
    const double* const delaunay_args) noexcept {
  /* nullify geopotential coeffs */
  mcs.clear();

  /* compute GMST using IAU 2006/2000A [rad] */
  const double gmst = dso::gmst(mjdtt, mjdut1);

  /* compute six-vector of multipliers ni from Delaunay vars */
  double __dargs[6];
  const double* __restrict__ f = dso::delaunay2doodson(delaunay_args, gmst, __dargs);

  /* iterate through individual constituents */
  for (const auto& wave : atlas().waves()) {
    /* compute angle: θ(f) = Σ(i=1,6) n(i)*β(i) */
    const double arg = wave.wave().doodson().argument(f);
    const double carg = std::cos(arg);
    const double sarg = std::sin(arg);
    mcs.Cnm() += wave.stokes_cos().Cnm() * carg + wave.stokes_sin().Cnm() * sarg;
    mcs.Snm() += wave.stokes_cos().Snm() * carg + wave.stokes_sin().Snm() * sarg;
  }

  return 0;
}

//int dso::OceanTide::stokes_coeffs(const dso::MjdEpoch& mjdtt,
//    const dso::MjdEpoch& mjdut1,
//    const double* const delaunay_args) noexcept {
//  /* nullify geopotential coeffs */
//  mcs.clear();
//
//  /* compute GMST using IAU 2006/2000A [rad] */
//  const double gmst = dso::gmst(mjdtt, mjdut1);
//
//  /* compute six-vector of multipliers ni from Delaunay vars */
//  double __dargs[6];
//  dso::delaunay2doodson(delaunay_args, gmst, __dargs);
//  Eigen::Matrix<double, 6, 1> theta = Eigen::Map<Eigen::Matrix<double, 6, 1>>(__dargs);
//
//  /* compute theta angles for all waves (fx1) */
//  const auto thetaf = doodson_matrix() * thetaf;
//
//  /* factors for cos and sin (kx) */
//  const auto factorCos = admittance_matrix().array().colwise() * thetaf.array().cos();
//  const auto factorSin = admittance_matrix().array().colwise() * thetaf.array().sin();
//
//  return 0;
//}
