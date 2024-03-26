#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "gravity.hpp"
#include "icgemio.hpp"
#include "planets.hpp"

using namespace costg;
constexpr const double GM_Moon    = 0.49028010560e13;

int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(
        stderr,
        "Usage: %s [00orbit_icrf.txt] [03directTideMoon_icrf.txt] [de421.bsp] [naif*.tls]\n",
        argv[0]);
    return 1;
  }

  /* Load CSPICE/NAIF Kernels */
  dso::load_spice_kernel(argv[3]);
  dso::load_spice_kernel(argv[4]);

  /* read orbit from input file */
  const auto orbvec = parse_orbit(argv[1]);
  printf("#Note: read %d data sets from input file\n", (int)orbvec.size());

  /* read accleration from input file */
  const auto accvec = parse_acceleration(argv[2]);
  printf("#Note: read %d data sets from input file\n", (int)accvec.size());

  /* compare results epoch by epoch */
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 1> rtb;
  Eigen::Matrix<double, 3, 3> g;
  auto acc = accvec.begin();
  for (const auto &in : orbvec) {
    /* GPSTime to TT */
    const auto t = in.epoch.gps2tai().tai2tt();
    /* get Moon position in ICRF */
    if (dso::planet_pos(dso::Planet::MOON, t, rtb)) {
      fprintf(stderr, "ERROR Failed to compute Moon position!\n");
      return 2;
    }
    /* compute third body acceleration for given epoch/position */
    a = dso::point_mass_acceleration(in.xyz, rtb, GM_Moon, g);
    /* get COSTG result */
    if (acc->epoch != in.epoch) {
      fprintf(stderr, "ERROR Failed to match epochs in input files\n");
      return 1;
    }
    printf("%d %.9f %.15e %.15e %.15e %.15e %.15e %.15e\n", in.epoch.imjd(),
           in.epoch.seconds(), acc->axyz(0), acc->axyz(1), acc->axyz(2), a(0),
           a(1), a(2));
    ++acc;
  }
}
