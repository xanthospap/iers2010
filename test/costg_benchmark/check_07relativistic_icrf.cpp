#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "eop.hpp"
#include "iersconst.hpp"
#include "planets.hpp"
#include "relativity.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

using namespace costg;
// constexpr const int formatD3Plot = 0;
constexpr const double TOLERANCE = 1e-11; /* [m/sec**2] */

int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(stderr,
            "Usage: %s [00orbit_icrf.txt] [07relativistic_icrf.txt] "
            "[de421.bsp] [naif*.tls]\n",
            argv[0]);
    return 1;
  }

  /* Load CSPICE/NAIF Kernels */
  dso::load_spice_kernel(argv[3]);
  dso::load_spice_kernel(argv[4]);

  /* read orbit from input file */
  const auto orbvec = parse_orbit(argv[1]);

  /* read acceleration from input file */
  const auto accvec = parse_acceleration(argv[2]);

  /* spit out a title for plotting */
  //if (formatD3Plot) {
  //  printf("mjd,sec,refval,val,component\n");
  //} else {
  //  printf("#title Relativistic Corrections (ephermeris: %s)\n",
  //         basename(argv[3]));
  //}

  /* compare results epoch by epoch */
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 6, 1> rsun;
  auto acc = accvec.begin();
  for (const auto &in : orbvec) {
    /* GPSTime to TT */
    // const auto t = in.epoch.gps2tai().tai2tt();
    const auto t = in.epoch;

    /* get Sun state in ICRF */
    if (dso::planet_state(dso::Planet::SUN, t.gps2tai().tai2tt(), rsun)) {
      fprintf(stderr, "ERROR Failed to compute Sun position!\n");
      return 2;
    }

    /* compute relativistic acceleration [m/sec^2] */
    Eigen::Matrix<double, 6, 1> rsat;
    rsat.head<3>() = in.xyz;
    rsat.tail<3>() = in.vxyz;
    a = dso::iers2010_relativistic_acceleration(rsat, rsun);

    /* get COSTG result */
    if (acc->epoch != in.epoch) {
      fprintf(stderr, "ERROR Failed to match epochs in input files\n");
      return 1;
    }

    //if (formatD3Plot) {
    //  printf("%d,%.9f,%.17e,%.17e,X\n", in.epoch.imjd(),
    //         in.epoch.seconds().seconds(), acc->axyz(0), a(0));
    //  printf("%d,%.9f,%.17e,%.17e,Y\n", in.epoch.imjd(),
    //         in.epoch.seconds().seconds(), acc->axyz(1), a(1));
    //  printf("%d,%.9f,%.17e,%.17e,Z\n", in.epoch.imjd(),
    //         in.epoch.seconds().seconds(), acc->axyz(2), a(2));
    //} else {
    //  printf("%d %.9f %.17e %.17e %.17e %.17e %.17e %.17e\n", in.epoch.imjd(),
    //         in.epoch.seconds().seconds(), acc->axyz(0), acc->axyz(1),
    //         acc->axyz(2), a(0), a(1), a(2));
    //}

    assert(std::abs(acc->axyz(0) - a(0)) < TOLERANCE);
    assert(std::abs(acc->axyz(1) - a(1)) < TOLERANCE);
    assert(std::abs(acc->axyz(2) - a(2)) < TOLERANCE);

    ++acc;
  }

  return 0;
}
