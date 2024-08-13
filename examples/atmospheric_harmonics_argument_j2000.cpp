#include "atmospheric_tides.hpp"
#include "doodson.hpp"
#include "geodesy/units.hpp"
#include "iau.hpp"
#include "eop.hpp"
#include "datetime/calendar.hpp"
#include "fundarg.hpp"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage %s [eopc]\n", argv[0]);
    return 1;
  }

  /* read EOPS (we will need DUT1) */
  dso::EopSeries eop;
  {
    auto imjd = dso::MjdEpoch::j2000_mjd();
    const auto t1 = dso::MjdEpoch(imjd - 2);
    const auto t2 = dso::MjdEpoch(imjd + 3);
    /* parse EOP values to the EopSeries instance for any t >= t1 and t < t2 */
    if (dso::parse_iers_C04(argv[1], t1, t2, eop)) {
      fprintf(stderr, "ERROR Failed parsing eop file\n");
      return 1;
    }
  }

  auto t = dso::MjdEpoch::j2000_mjd();

  /* compute gmst using an approximate value for UT1 (linear interpolation) */
  double dut1_approx;
  eop.approx_dut1(t, dut1_approx);

  /* compute GMST using IAU 2006/2000A [rad] */
  auto tut1 = t.tt2ut1(dut1_approx);
  const double gmst = dso::gmst(t, tut1);

  /* compute fundamental (Delaunay) arguments fot t */
  double fargs[5];
  dso::fundarg(t, fargs);

  /* fundamental (Delaunay) to Doodson arguments */
  double dargs[6];
  dso::delaunay2doodson(fargs, gmst, dargs);

  char buf[16];
  for (const auto &harmonic : dso::TidalConstituentsArray) {
    /* compute angle: θ(f) = Σ(i=1,6) n(i)*β(i) + χ */
    const double arg =
        harmonic._d.argument(dargs) + harmonic._d.pifactor() * dso::DPI/2;
    printf("Harmonic %s argument at J2000 : %.3f [deg]\n",
           harmonic._d.str(buf, true), dso::rad2deg(arg));
  }

  return 0;
}
