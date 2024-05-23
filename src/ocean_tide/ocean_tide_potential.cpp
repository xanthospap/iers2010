#include "doodson.hpp"
#include "iau.hpp"
#include "ocean_tide.hpp"

int dso::OceanTide::stokes_coeffs(int max_degree, int max_order,
                                  const dso::MjdEpoch &mjdtt,
                                  const dso::MjdEpoch &mjdut1,
                                  const double *const delaunay_args) noexcept {
  /* check d/o */
  assert(max_degree >= max_order);

  /* nullify geopotential coeffs */
  mcs.clear();

  /* compute GMST using IAU 2006/2000A [rad] */
  const double gmst = dso::gmst(mjdtt, mjdut1);

  /* compute six-vector of multipliers ni from Delaunay vars */
  double __dargs[6];
  const double *__restrict__ f =
      dso::delaunay2doodson(delaunay_args, gmst, __dargs);

  /* iterate through individual constituents */
  for (const auto &wave : mwaves) {
    /* compute angle: θ(f) = Σ(i=1,6) n(i)*β(i) */
    const double arg = wave.doodson().argument(f);
    const double carg = std::cos(arg);
    const double sarg = std::sin(arg);
    mcs.Cnm() +=
        wave.cos_coeffs().Cnm() * carg + wave.sin_coeffs().Cnm() * sarg;
    mcs.Snm() +=
        wave.cos_coeffs().Snm() * carg + wave.sin_coeffs().Snm() * sarg;
  }

  return 0;
}
