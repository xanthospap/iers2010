#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "earth_rotation.hpp"
#include "eigen3/Eigen/Eigen"
#include "eop.hpp"
#include "fundarg.hpp"
#include "gravity.hpp"
#include "iau.hpp"
#include "icgemio.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "sofa.h"
#include <cassert>

constexpr const double PTOLERANCE = 1e-13; /* [m/sec**2] */
constexpr const double VTOLERANCE = 1e-13; /* [m/sec**2] */
using namespace costg;

/* Celestial-to-terrestrial via SOFA, i.e. [TRS] = rc2t * [CRS] */
Eigen::Matrix<double, 3, 3> sofa(double tt1, double tt2, double ut1, double ut2,
                                 double xp, double yp) {
  double rc2t[3][3];
  iauC2t06a(tt1, tt2, ut1, ut2, xp, yp, rc2t);
  /* copy to output matrix */
  Eigen::Matrix<double, 3, 3> R;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      R(i, j) = rc2t[i][j];
    }
  }
  return R;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr,
            "Usage: %s [00orbit_itrf.txt] [00orbit_icrf.txt] "
            "[eopc04_14_IAU2000.62-now]\n",
            argv[0]);
    return 1;
  }

  /* read orbit from input file */
  const auto torbvec = parse_orbit(argv[1]);

  /* read orbit from input file */
  const auto corbvec = parse_orbit(argv[2]);

  /* create an instance to hold EOP series */
  dso::EopSeries eops;

  /* parse EOP values to the EopSeries instance for any t >= t1 and t < t2 */
  if (dso::parse_iers_C04(
          argv[3],
          torbvec[0].epoch.add_seconds(dso::FractionalSeconds(-2 * 86400e0)),
          torbvec[0].epoch.add_seconds(dso::FractionalSeconds(2 * 86400e0)),
          eops)) {
    fprintf(stderr, "ERROR Failed parsing eop file\n");
    return 1;
  }

  /* regularize the series (UT1 and LOD) */
  eops.regularize();

  /* compare results epoch by epoch */
  Eigen::Matrix<double, 3, 3> dRdt;
  auto tin = torbvec.begin();
  for (const auto &cin : corbvec) {
    /* get current time in TT */
    const auto tt = cin.epoch.gps2tt();

    /* interpolate EOPS */
    dso::EopRecord eopr;
    if (dso::EopSeries::out_of_bounds(eops.interpolate(tt, eopr))) {
      fprintf(stderr, "Failed to interpolate: Epoch is out of bounds!\n");
      return 1;
    }

    /* compute gmst using an approximate value for UT1 (linear interpolation) */
    double dut1_approx;
    eops.approx_dut1(tt, dut1_approx);
    const double gmst = dso::gmst(tt, tt.tt2ut1(dut1_approx));

    /* compute fundamental arguments at given epoch */
    double fargs[14];
    dso::fundarg(tt, fargs);

    /* add libration effect [micro as] */
    {
      double dxp, dyp, dut1, dlod;
      dso::deop_libration(fargs, gmst, dxp, dyp, dut1, dlod);
      eopr.xp() += dxp * 1e-6;
      eopr.yp() += dyp * 1e-6;
      eopr.dut() += dut1 * 1e-6;
      eopr.lod() += dlod * 1e-6;
    }

    /* add ocean tidal effect [micro as] */
    {
      double dxp, dyp, dut1, dlod;
      dso::deop_ocean_tide(fargs, gmst, dxp, dyp, dut1, dlod);
      eopr.xp() += dxp * 1e-6;
      eopr.yp() += dyp * 1e-6;
      eopr.dut() += dut1 * 1e-6;
      eopr.lod() += dlod * 1e-6;
    }

    /* get the rotation matrix, i.e. [TRS] = q * [CRS] */
    Eigen::Quaterniond q = dso::c2i06a(tt, eopr, dRdt);

    // Eigen::Matrix<double, 3, 3> Rs;
    //{
    //   const double jd1_tt = tt.imjd() + dso::MJD0_JD;
    //   const double jd2_tt = tt.fractional_days().days();
    //   const auto ut = tt.tt2ut1(eopr.dut());
    //   const double jd1_ut = ut.imjd() + dso::MJD0_JD;
    //   const double jd2_ut = ut.fractional_days().days();
    //   Rs = sofa(jd1_tt, jd2_tt, jd1_ut, jd2_ut, dso::sec2rad(eopr.xp()),
    //             dso::sec2rad(eopr.yp()));
    // }

    /* transform position */
    const Eigen::Vector3d mtpos = q * cin.xyz;
    /* transform velocity */
    const auto mtvel = q * cin.vxyz + dRdt * cin.xyz;
    /* transform acceleration */
    const auto mtacc = q * cin.axyz;

    {
      printf("R=\n");
      const auto Q = q.toRotationMatrix();
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          printf("+%.3e ", Q(i, j));
        }
        printf("\n");
      }
      printf("dRdt=\n");
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          printf("+%.3e ", dRdt(i, j));
        }
        printf("\n");
      }
    }

    /* get COSTG result */
    if (tin->epoch != cin.epoch) {
      fprintf(stderr, "ERROR Faile to match epochs in input files\n");
      return 1;
    }

    printf("dP = %+.6f %+.6f %+.6f\n", mtpos(0) - tin->xyz(0),
           mtpos(1) - tin->xyz(1), mtpos(2) - tin->xyz(2));
    printf("dV = %+.6f %+.6f %+.6f\n", mtvel(0) - tin->vxyz(0),
           mtvel(1) - tin->vxyz(1), mtvel(2) - tin->vxyz(2));
    printf("dA = %+.6f %+.6f %+.6f\n", mtacc(0) - tin->axyz(0),
           mtacc(1) - tin->axyz(1), mtacc(2) - tin->axyz(2));

    // assert(std::abs(acc->axyz(0) - a(0)) < TOLERANCE);
    // assert(std::abs(acc->axyz(1) - a(1)) < TOLERANCE);
    // assert(std::abs(acc->axyz(2) - a(2)) < TOLERANCE);

    ++tin;
  }

  return 0;
}
