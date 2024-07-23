#include "fundarg.hpp"
#include "solid_earth_tide.hpp"
#include "earth_rotation.hpp"
#include "planets.hpp"

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "Usage: [eopc04.1962-now] [de421.bsp] [naif*.tls]\n");
    return 1;
  }

  /* Approx ITRF coordinates for site DIOB */
  Eigen::Matrix<double, 3, 1> rsta;
  rsta << 4595212.468e0, 2039473.691e0, 3912617.891e0;

  /* start and end epochs in TT */
  MjdEpoch start(year(2023), month(7), day_of_month(19));
  MjdEpoch end(year(2023), month(7), day_of_month(29));

  /* create an instance to hold EOP series */
  EopSeries eops;

  /* parse EOP values to the EopSeries instance for any t >= t1 and t < t2 */
  if (parse_iers_C04(argv[1], start - modified_julian_day(2),
                     end + modified_julian_day(2), eops)) {
    fprintf(stderr, "ERROR Failed parsing eop file\n");
    return 1;
  }

  /* Load CSPICE/NAIF Kernels */
  load_spice_kernel(argv[2]);
  load_spice_kernel(argv[3]);

  /* setup a SolidEarthTide instance */
  SolidEarthTide set(3.986004418e14, 6378136.6e0, 4.9048695e12,
                     1.32712442099e20);
  
  Eigen::Matrix<double, 3, 1> rmon, rsun;
  auto tt = start;
  while (tt<end) {
    
    /* Fundamental arguments (IERS 2003) */
    const double fa[14] = {
        fal03(tt),  falp03(tt), faf03(tt),  fad03(tt),  faom03(tt),
        fame03(tt), fave03(tt), fae03(tt),  fama03(tt), faju03(tt),
        fasa03(tt), faur03(tt), fane03(tt), fapa03(tt),
    };

    /* EOPS at time of request */
    EopRecord eop;
    if (EopSeries::out_of_bounds(eops.interpolate(tt, eop))) {
      fprintf(stderr, "Failed to interpolate: Epoch is out of bounds!\n");
      return 1;
    }

    /* get Sun+Moon position in ICRF */
    if (planet_pos(Planet::SUN, tt, rsun)) {
      fprintf(stderr, "ERROR Failed to compute Sun position!\n");
      return 2;
    }
    if (planet_pos(Planet::MOON, tt, rmon)) {
      fprintf(stderr, "ERROR Failed to compute Sun position!\n");
      return 2;
    }

    /* transformation quaternion */
    const auto q = dso::itrs2gcrs_quaternion(tt, eop);
    rsun = q.inverse() * rsun;
    rmon = q.inverse() * rmon;

    /* compute deformation (cartesian) */
    auto dr = set.displacement(tt, tt.tt2ut1(eop.dut()), rsta, rmon, rsun, fa);

    printf("%.9f %.4f %.4f %.4f %.4f\n", tt.imjd() + tt.fractional_days(),
           dr(0), dr(1), dr(2), dr.norm());

    /* next date */
    tt.add_seconds(FractionalSeconds(30));
  }

  return 0;
}
