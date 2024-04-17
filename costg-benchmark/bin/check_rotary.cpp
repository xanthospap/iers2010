#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "gravity.hpp"
#include "icgemio.hpp"
#include "planets.hpp"
#include "iersconst.hpp"
#include "solid_earth_tide.hpp"
#include "eop.hpp"

using namespace costg;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(
        stderr,
        "Usage: %s [00orbit_icrf.txt] [00orbit_itrf.txt] [01earthRotation_rotaryMatrix.txt]\n",
        argv[0]);
    return 1;
  }

  /* read orbit from input file */
  const auto orbvec_cel = parse_orbit(argv[1]);
  
  /* read orbit from input file */
  const auto orbvec_ter = parse_orbit(argv[2]);
  
  /* read rotary matrix (GCRS to ITRS) from input file */
  const auto rotvec = parse_rotary(argv[3]);

  Eigen::Matrix<double, 3, 1> r;
  auto ritrf = orbvec_ter.begin();
  auto rot = rotvec.begin();
  for (const auto &in : orbvec_cel) {

    /* rotation matrix */
    const Eigen::Matrix<double,3,3> R = rot->R;
    assert(rot->epoch == in.epoch);
    
    /* GCRF to ITRF */
    r = R * in.xyz;
    
    /* get COSTG result */
    if (ritrf->epoch != in.epoch) {
      fprintf(stderr, "ERROR Failed to match epochs in input files\n");
      return 1;
    }

    printf("%d %.9f %.15e %.15e %.15e %.15e %.15e %.15e\n", in.epoch.imjd(),
           in.epoch.seconds(), ritrf->xyz(0), ritrf->xyz(1), ritrf->xyz(2), r(0),
           r(1), r(2));
    ++rot;
    ++ritrf;
  }
}
