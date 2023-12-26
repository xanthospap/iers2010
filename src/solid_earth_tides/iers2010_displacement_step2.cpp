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

constexpr const std::array<Table7, 16> Tables7 = {
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
#ifdef DEBUG
  //printf("T=%19.12f FHR=%15.10f/%15.10f\n", 
  //mjdtt.jcenturies_sinceJ2000(), mjdtt.seconds()/3600e0, mjdut1.seconds()/3600e0);
  //printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", rad2deg(f[0]),
  //       rad2deg(f[1]), rad2deg(f[2]), rad2deg(f[3]), rad2deg(f[4]),
  //       rad2deg(f[5]));
#endif

  /* displacement vector (to be returned) */
  Eigen::Matrix<double, 3, 1> dr = Eigen::Matrix<double, 3, 1>::Zero();
#ifdef DEBUG
  Eigen::Matrix<double, 3, 1> ddr = Eigen::Matrix<double, 3, 1>::Zero();
  Eigen::Matrix<double, 3, 1> ldr = Eigen::Matrix<double, 3, 1>::Zero();
#endif

  /* Contributions from the diurnal band, for n=2 */
  for (const auto &tc : Tables7) {
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
#ifdef DEBUG
      //printf("[dirnal] theta: %10.6f SIN=  %15.9f  COS=%15.9f\n", arg, std::sin(arg), std::cos(arg));
      ddr(2) += (__sa * tc.dRip + __ca * tc.dRop) * __s2phi;
      ddr(0) += (__ca * tc.dTip - __sa * tc.dTop) * __sphi;
      ddr(1) += (__sa * tc.dTip + __ca * tc.dTop) * __c2phi;
      //printf("Correction     : %15.9f %15.9f %15.9f\n", ddr(0), ddr(1), ddr(2));
#endif
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
#ifdef DEBUG
      ldr(2) += (1.50e0 * __sphi * __sphi - .5e0) * (tc.dRip * __ca + tc.dRop * __sa);
               
      ldr(0) += 0e0;
      ldr(1) += (tc.dTip * __ca + tc.dTop * __sa) * __s2phi;
#endif
    }
  }

  /* [mm] to [m] */
  dr *= 1e-3;
#ifdef DEBUG
  const auto R = dso::geodetic2lvlh(tsta.msph.lat(), tsta.msph.lon());
  ddr = R * ddr;
  ldr = R * ldr;
  ddr *= 1e-3;
  ldr *= 1e-3;
  printf("STEP2 Diurnal : %10.6f %10.6f %10.6f\n", ddr(0), ddr(1), ddr(2));
  printf("STEP2 Long    : %10.6f %10.6f %10.6f\n", ldr(0), ldr(1), ldr(2));
  //printf("STEP2 long    : %10.6f %10.6f %10.6f\n", ldr(0), ldr(1), ldr(2));
#endif
  return dr;
}
