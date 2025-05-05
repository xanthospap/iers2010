#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "earth_rotation.hpp"
#include "eigen3/Eigen/Eigen"
#include "eop.hpp"
#include "fundarg.hpp"
#include "gravity.hpp"
#include "iau.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "sofa.h"
#include <cassert>

constexpr const double STOLERANCE = 1e-1; /* in [arc "] */
constexpr const double TTOLERANCE = 1e+1; /* in ["] */

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr,
            "Usage: %s [01earthRotation_interpolatedEOP.txt] "
            "[eopc04_14_IAU2000.62-now]\n",
            argv[0]);
    return 1;
  }

  /* read EOPS from input file */
  const auto ceops = costg::parse_eops(argv[1]);

  /* create an instance to hold EOP series */
  dso::EopSeries eops;

  /* parse EOP values to the EopSeries instance for any t >= t1 and t < t2 */
  if (dso::parse_iers_C04(
          argv[2],
          ceops[0].epoch.add_seconds(dso::FractionalSeconds(-2 * 86400e0)),
          ceops[0].epoch.add_seconds(dso::FractionalSeconds(2 * 86400e0)),
          eops)) {
    fprintf(stderr, "ERROR Failed parsing eop file\n");
    return 1;
  }

  /* regularize the series (UT1 and LOD) */
  eops.regularize();

  /* check */
  dso::EopRecord meop;
  double fargs[14];
  double Xcip, Ycip;

  for (const auto &ceop : ceops) {
    /* interpolate EOPS */
    if (dso::EopSeries::out_of_bounds(
            eops.interpolate(ceop.epoch.gps2tt(), meop))) {
      fprintf(stderr, "Failed to interpolate: Epoch is out of bounds!\n");
      return 1;
    }

    /* compute (X,Y) CIP and fundamental arguments */
    dso::xycip06a(ceop.epoch.gps2tt(), Xcip, Ycip, fargs);

    /* use fundamental arguments to compute s */
    const double s = dso::s06(ceop.epoch.gps2tt(), Xcip, Ycip, fargs);

    /* apply CIP corrections */
    Xcip += dso::sec2rad(meop.dX());
    Ycip += dso::sec2rad(meop.dY());

    /* compute gmst using an approximate value for UT1 (linear interpolation) */
    double dut1_approx;
    eops.approx_dut1(ceop.epoch.gps2tt(), dut1_approx);
    const double gmst =
        dso::gmst(ceop.epoch.gps2tt(), ceop.epoch.gps2tt().tt2ut1(dut1_approx));

    /* add libration effect [micro as] */
    {
      double dxp, dyp, dut1, dlod;
      dso::deop_libration(fargs, gmst, dxp, dyp, dut1, dlod);
      meop.xp() += dxp * 1e-6;
      meop.yp() += dyp * 1e-6;
      meop.dut() += dut1 * 1e-6;
      meop.lod() += dlod * 1e-6;
    }

    /* add ocean tidal effect [micro as] */
    {
      double dxp, dyp, dut1, dlod;
      dso::deop_ocean_tide(fargs, gmst, dxp, dyp, dut1, dlod);
      meop.xp() += dxp * 1e-6;
      meop.yp() += dyp * 1e-6;
      meop.dut() += dut1 * 1e-6;
      meop.lod() += dlod * 1e-6;
    }

    /* de-regularize */
    {
      double ut1_cor;
      double lod_cor;
      double omega_cor;

      dso::deop_zonal_tide(fargs, ut1_cor, lod_cor, omega_cor);
      /* apply (note: microseconds to seconds) */
      meop.dut() += (ut1_cor * 1e-6);
      meop.lod() += (lod_cor * 1e-6);
    }

    /* check 1. diffs in ["] "*/
    // printf("%d %.9f %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
    // ceop.epoch.imjd(), ceop.epoch.seconds().seconds(),
    // dso::rad2sec(ceop.xp - dso::sec2rad(meop.xp())),
    // dso::rad2sec(ceop.yp - dso::sec2rad(meop.yp())),
    // dso::rad2sec(ceop.sp - sp00(ceop.epoch.gps2tt())),
    // ceop.dUT1 - meop.dut(), ceop.LOD - meop.lod(),
    // dso::rad2sec(ceop.X - Xcip), dso::rad2sec(ceop.Y - Ycip),
    // dso::rad2sec(ceop.s - s));

    /* check 2. this lib interpolation results, same units as COSTG benchmark */
    printf("%d %.9f %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
           ceop.epoch.imjd(), ceop.epoch.seconds().seconds(),
           dso::sec2rad(meop.xp()), dso::sec2rad(meop.yp()),
           sp00(ceop.epoch.gps2tt()), meop.dut(), meop.lod(), Xcip, Ycip, s);

    /* assert diffs */
    assert(dso::rad2sec(std::abs(ceop.xp - dso::sec2rad(meop.xp()))) <
           STOLERANCE);
    assert(dso::rad2sec(std::abs(ceop.yp - dso::sec2rad(meop.yp()))) <
           STOLERANCE);
    assert(dso::rad2sec(std::abs(ceop.sp - sp00(ceop.epoch.gps2tt()))) <
           STOLERANCE);
    assert(std::abs(ceop.dUT1 - meop.dut()) < TTOLERANCE);
    assert(std::abs(ceop.LOD - meop.lod()) < TTOLERANCE);
    assert(dso::rad2sec(std::abs(ceop.X - Xcip)) < STOLERANCE);
    assert(dso::rad2sec(std::abs(ceop.Y - Ycip)) < STOLERANCE);
    assert(dso::rad2sec(std::abs(ceop.s - s)) < STOLERANCE);
  }

  return 0;
}