#include "eop.hpp"
#include "fundarg.hpp" /* for fundarg */
#include "iau.hpp"     /* for gmst */
#include <cassert>
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage %s [EOP C04 FILE]\n", argv[0]);
    return 1;
  }

  /* 28 December 2016 */
  dso::MjdEpoch t1(57750);
  /* 6 January 2017 (not inclusive) */
  dso::MjdEpoch t2(57759);

  /* create an instance to hold EOP series */
  dso::EopSeries eop;

  /* parse EOP values to the EopSeries instance for any t >= t1 and t < t2 */
  if (dso::parse_iers_C04(argv[1], t1, t2, eop)) {
    fprintf(stderr, "ERROR Failed parsing eop file\n");
    return 1;
  }

  /* check entries parsed */
  assert(eop.num_entries() == 9);

  /* regularize the series (UT1 and LOD) */
  eop.regularize();

  /* storage for fundamental arguments */
  double fargs[5];

  /* interpolate every 30 seconds; hold results in a EopRecord instance */
  dso::EopRecord eops;
  auto t = eop.first_epoch();
  while (t < eop.last_epoch()) {
    /* interpolate */
    if (dso::EopSeries::out_of_bounds(eop.interpolate(t, eops))) {
      fprintf(stderr, "Failed to interpolate: Epoch is out of bounds!\n");
      return 1;
    }

    /* print raw interpolation results */
    printf("[RAW] %.9f %.12f %.12f %.12f %.12f\n", eops.t().as_mjd(), eops.xp(),
           eops.yp(), eops.dut(), eops.lod());

    /* compute gmst using an approximate value for UT1 (linear interpolation) */
    double dut1_approx;
    eop.approx_dut1(t, dut1_approx);
    const double gmst = dso::gmst(t, t.tt2ut1(dut1_approx));

    /* compute fundamental arguments at given epoch */
    dso::fundarg(t, fargs);

    /* add libration effect [micro as] */
    {
      double dxp, dyp, dut1, dlod;
      dso::deop_libration(fargs, gmst, dxp, dyp, dut1, dlod);
      eops.xp() += dxp * 1e-6;
      eops.yp() += dyp * 1e-6;
      eops.dut() += dut1 * 1e-6;
      eops.lod() += dlod * 1e-6;
    }

    /* add ocean tidal effect [micro as] */
    {
      double dxp, dyp, dut1, dlod;
      dso::deop_ocean_tide(fargs, gmst, dxp, dyp, dut1, dlod);
      eops.xp() += dxp * 1e-6;
      eops.yp() += dyp * 1e-6;
      eops.dut() += dut1 * 1e-6;
      eops.lod() += dlod * 1e-6;
    }

    /* print corrected interpolation results */
    printf("[COR] %.9f %.12f %.12f %.12f %.12f\n", eops.t().as_mjd(), eops.xp(),
           eops.yp(), eops.dut(), eops.lod());
    t.add_seconds(dso::FractionalSeconds(30));
  }

  /* optional: de-regularize eop series */
  eop.deregularize();

  return 0;
}
