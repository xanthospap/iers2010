#include "solid_earth_tide.hpp"
#include "doodson.hpp"
#include "iau.hpp"
#include <array>
#include <cmath>

namespace {
/** A struct to hold lines of Tables 7a and 7b in the IERS2010 standards */
struct Table7 {
  /** Doodson number */
  dso::DoodsonConstituent md;
  /** ΔR in-phase for this doodson number [mm] */
  double dRip;
  /** ΔR out-of-phase for this doodson number [mm] */
  double dRop;
  /** ΔT in-phase for this doodson number [mm] */
  double dTip;
  /** ΔT out-of-phase for this doodson number [mm] */
  double dTop;
}; /* Table7 */

/** @brief Diurnal and long waves for computing solid earth todes Step-2 
 * corrections.
 * 
 * This is a translation of the Tables 7a and 7b of the IERS2010 standards.
 */
constexpr const std::array<Table7, 16> Tables7 = {
    /* Diurnal tidal waves */
    {{/*135655*/ {1, -2, 0, 1, 0, 0}, -8.0000e-02, +0.0000e+00, -1.0000e-02,
      +1.0000e-02},
     {/*145545*/ {1, -1, 0, 0, -1, 0}, -1.0000e-01, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*145555*/ {1, -1, 0, 0, 0, 0}, -5.1000e-01, +0.0000e+00, -2.0000e-02,
      +3.0000e-02},
     {/*155655*/ {1, 0, 0, 1, 0, 0}, +6.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*162556*/ {1, 1, -3, 0, 0, 1}, -6.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*163555*/ {1, 1, -2, 0, 0, 0}, -1.2300e+00, -7.0000e-02, +6.0000e-02,
      +1.0000e-02},
     {/*165545*/ {1, 1, 0, 0, -1, 0}, -2.2000e-01, +1.0000e-02, +1.0000e-02,
      +0.0000e+00},
     {/*165555*/ {1, 1, 0, 0, 0, 0}, +1.2000e+01, -7.8000e-01, -6.7000e-01,
      -3.0000e-02},
     {/*165565*/ {1, 1, 0, 0, 1, 0}, +1.7300e+00, -1.2000e-01, -1.0000e-01,
      +0.0000e+00},
     {/*166554*/ {1, 1, 1, 0, 0, -1}, -5.0000e-01, -1.0000e-02, +3.0000e-02,
      +0.0000e+00},
     {/*167555*/ {1, 1, 2, 0, 0, 0}, -1.1000e-01, +1.0000e-02, +1.0000e-02,
      +0.0000e+00},
      /* Long-period tidal waves */
     {/*55565*/ {0, 0, 0, 0, 1, 0}, +4.7000e-01, +1.6000e-01, +2.3000e-01,
      +7.0000e-02},
     {/*57555*/ {0, 0, 2, 0, 0, 0}, -2.0000e-01, -1.1000e-01, -1.2000e-01,
      -5.0000e-02},
     {/*65455*/ {0, 1, 0, -1, 0, 0}, -1.1000e-01, -9.0000e-02, -8.0000e-02,
      -4.0000e-02},
     {/*75555*/ {0, 2, 0, 0, 0, 0}, -1.3000e-01, -1.5000e-01, -1.1000e-01,
      -7.0000e-02},
     {/*75565*/ {0, 2, 0, 0, 1, 0}, -5.0000e-02, -6.0000e-02, -5.0000e-02,
      -3.0000e-02}}};

/** @brief Diurnal and long waves for computing solid earth todes Step-2 
 * corrections.
 * 
 * The coefficients and waves here are 31 for the diurnal band and not 
 * 11 as in the IERS2010 standards. This table is extracted from the 
 * DEAHNTTIDEINEL.F subroutine, published by IERS.
 */
constexpr const std::array<Table7, 36> Tables7Extended = {
    /* Diurnal tidal waves */
    {{/*125755*/ {1, -3, 0, 2, 0, 0}, -1.0000e-02, -1.0000e-02, +0.0000e+00,
      +0.0000e+00},
     {/*127555*/ {1, -3, 2, 0, 0, 0}, -1.0000e-02, -1.0000e-02, +0.0000e+00,
      +0.0000e+00},
     {/*135645*/ {1, -2, 0, 1, -1, 0}, -2.0000e-02, -1.0000e-02, +0.0000e+00,
      +0.0000e+00},
     {/*135655*/ {1, -2, 0, 1, 0, 0}, -8.0000e-02, +0.0000e+00, +1.0000e-02,
      +1.0000e-02},
     {/*137455*/ {1, -2, 2, -1, 0, 0}, -2.0000e-02, -1.0000e-02, +0.0000e+00,
      +0.0000e+00},
     {/*145545*/ {1, -1, 0, 0, -1, 0}, -1.0000e-01, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*145555*/ {1, -1, 0, 0, 0, 0}, -5.1000e-01, +0.0000e+00, -2.0000e-02,
      +3.0000e-02},
     {/*147555*/ {1, -1, 2, 0, 0, 0}, +1.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*153655*/ {1, 0, -2, 1, 0, 0}, +1.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*155455*/ {1, 0, 0, -1, 0, 0}, +2.0000e-02, +1.0000e-02, +0.0000e+00,
      +0.0000e+00},
     {/*155655*/ {1, 0, 0, 1, 0, 0}, +6.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*155665*/ {1, 0, 0, 1, 1, 0}, +1.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*157455*/ {1, 0, 2, -1, 0, 0}, +1.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*162556*/ {1, 1, -3, 0, 0, 1}, -6.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*163565*/ {1, 1, -2, 0, 1, 0}, +1.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*163555*/ {1, 1, -2, 0, 0, 0}, -1.2300e+00, -7.0000e-02, +6.0000e-02,
      +1.0000e-02},
     {/*164554*/ {1, 1, -1, 0, 0, -1}, +2.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*164556*/ {1, 1, -1, 0, 0, 1}, +4.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*165545*/ {1, 1, 0, 0, -1, 0}, -2.2000e-01, +1.0000e-02, +1.0000e-02,
      +0.0000e+00},
     {/*165555*/ {1, 1, 0, 0, 0, 0}, +1.2000e+01, -7.8000e-01, -6.7000e-01,
      -3.0000e-02},
     {/*165565*/ {1, 1, 0, 0, 1, 0}, +1.7300e+00, -1.2000e-01, -1.0000e-01,
      +0.0000e+00},
     {/*165575*/ {1, 1, 0, 0, 2, 0}, -4.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*166554*/ {1, 1, 1, 0, 0, -1}, -5.0000e-01, -1.0000e-02, +3.0000e-02,
      +0.0000e+00},
     {/*166556*/ {1, 1, 1, 0, 0, 1}, +1.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*166564*/ {1, 1, 1, 0, 1, -1}, -1.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*167355*/ {1, 1, 2, -2, 0, 0}, -1.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*167555*/ {1, 1, 2, 0, 0, 0}, -1.1000e-01, +1.0000e-02, +1.0000e-02,
      +0.0000e+00},
     {/*173655*/ {1, 2, -2, 1, 0, 0}, -1.0000e-02, +0.0000e+00, +0.0000e+00,
      +0.0000e+00},
     {/*175455*/ {1, 2, 0, -1, 0, 0}, -2.0000e-02, +2.0000e-02, +0.0000e+00,
      +1.0000e-02},
     {/*185555*/ {1, 3, 0, 0, 0, 0}, +0.0000e+00, +1.0000e-02, +0.0000e+00,
      +1.0000e-02},
     {/*185565*/ {1, 3, 0, 0, 1, 0}, +0.0000e+00, +1.0000e-02, +0.0000e+00,
      +0.0000e+00},
     /* Long-period tidal waves */
     {/*55565*/ {0, 0, 0, 0, 1, 0}, +4.7000e-01, +1.6000e-01, +2.3000e-01,
      +7.0000e-02},
     {/*57555*/ {0, 0, 2, 0, 0, 0}, -2.0000e-01, -1.1000e-01, -1.2000e-01,
      -5.0000e-02},
     {/*65455*/ {0, 1, 0, -1, 0, 0}, -1.1000e-01, -9.0000e-02, -8.0000e-02,
      -4.0000e-02},
     {/*75555*/ {0, 2, 0, 0, 0, 0}, -1.3000e-01, -1.5000e-01, -1.1000e-01,
      -7.0000e-02},
     {/*75565*/ {0, 2, 0, 0, 1, 0}, -5.0000e-02, -6.0000e-02, -5.0000e-02,
      -3.0000e-02}}};
} /* unnamed namespace */

Eigen::Matrix<double, 3, 1> dso::SolidEarthTide::step2_displacement(
    const dso::MjdEpoch &mjdtt, const dso::MjdEpoch &mjdut1,
    [[maybe_unused]]const Eigen::Matrix<double, 3, 1> &rsta,
    const dso::SolidEarthTide::PointSphericalTrigs &tsta,
    const double *const delaunay_args) noexcept {

  /* trigs */
  const double __sphi = tsta.mslat;
  const double __s2phi = 2e0 * tsta.mslat * tsta.mclat;
  const double __c2phi = tsta.mclat * tsta.mclat - tsta.mslat * tsta.mslat;
  const double lambda = tsta.msph.lon();

  /* compute GMST using IAU 2006/2000A [rad] */
  const double gmst = dso::gmst(mjdtt, mjdut1);
  
  /* get Doodson arguments from Delaunay args. */
  double __dargs[6];
  const double *__restrict__ f =
      dso::delaunay2doodson(delaunay_args, gmst, __dargs);

  /* displacement vector (to be returned) */
  Eigen::Matrix<double, 3, 1> dr = Eigen::Matrix<double, 3, 1>::Zero();

  /* Contributions from the diurnal band, for n=2 */
  for (const auto &tc : Tables7Extended) {
    /* θf */
    const double arg = tc.md.argument(f);
    if (tc.md.operator()(0) == 1) {
      /* diurnal band ang = (θf + λ), Eqs. 12 */
      const double ang = dso::anp(arg + lambda);
      const double __sa = std::sin(ang);
      const double __ca = std::cos(ang);
      /* east */
      dr(0) += (__ca * tc.dTip - __sa * tc.dTop) * __sphi;
      /* north */
      dr(1) += (__sa * tc.dTip + __ca * tc.dTop) * __c2phi;
      /* radial */
      dr(2) += (__sa * tc.dRip + __ca * tc.dRop) * __s2phi;
    } else {
      /* long-period band, Eqs. 13 */
      const double __sa = std::sin(arg);
      const double __ca = std::cos(arg);
      /* east */
      dr(0) += 0e0;
      /* north */
      dr(1) += (tc.dTip * __ca + tc.dTop * __sa) * __s2phi;
      /* radial */
      dr(2) += (1.50e0 * __sphi * __sphi - .5e0) *
               (tc.dRip * __ca + tc.dRop * __sa);
    }
  }

  /* [mm] to [m] */
  dr *= 1e-3;
  
  return dr;
}
