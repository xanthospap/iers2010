#include "atmospheric_tides.hpp"
#include "doodson.hpp"
#include "iau.hpp"

int dso::AtmosphericTides::stokes_coeffs(
    const dso::MjdEpoch &mjdtt, const dso::MjdEpoch &mjdut1,
    const double *const delaunay_args) noexcept {
  /* nullify geopotential coeffs */
  mcs.clear();
  
  /* compute GMST using IAU 2006/2000A [rad] */
  const double gmst = dso::gmst(mjdtt, mjdut1);

  /* compute six-vector of multipliers ni from Delaunay vars */
  double __dargs[6];
  const double *__restrict__ f =
      dso::delaunay2doodson(delaunay_args, gmst, __dargs);

  /* iterate through individual constituents */
  for (const auto &wave: mwaves) {
    /* compute angle: θ(f) = Σ(i=1,6) n(i)*β(i) */
    const double arg = wave.mdentry._d.argument(f) + wave.mdentry._d.pifactor();
    const double carg = std::cos(arg);
    const double sarg = std::sin(arg);
    mcs.Cnm() += wave.mCosCs.Cnm() * carg + wave.mSinCs.Cnm() * sarg;
    mcs.Snm() += wave.mCosCs.Snm() * carg + wave.mSinCs.Snm() * sarg;
  }

  return 0;
}
