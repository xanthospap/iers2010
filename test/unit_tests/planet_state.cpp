#include "planets.hpp"
#include <cstdio>
#include <datetime/dtfund.hpp>

int main(int argc, char *argv[]) {

  if (argc != 3) {
    fprintf(stderr, "Usage: %s [JPL EPHEMERIS <.bsp>] [JPL LEASP SEC FILE <.tls>]\n", argv[0]);
    return 1;
  }

  /* Load CSPICE/NAIF Kernels */
  dso::load_spice_kernel(argv[1]);
  dso::load_spice_kernel(argv[2]);

  constexpr const int numtests = 1000;

  /* For a series of random dates in the interval 01/01/1990 - 01/01/2030, 
   * get SUN/MOON pos and state (in GCRF), and compare position vectors. 
   * They should be the same.
   */
  Eigen::Matrix<double,3,1> pos;
  Eigen::Matrix<double,6,1> state;
  for (int i=0; i<numtests; i++) {
    dso::MjdEpoch t = dso::MjdEpoch::random(dso::modified_julian_day(47892),
                                  dso::modified_julian_day(62502));

    if (dso::planet_pos(dso::Planet::MOON, t, pos)) {
      fprintf(stderr, "ERROR. Failed getting Moon's position\n");
      return 5;
    }

    if (dso::planet_state(dso::Planet::MOON, t, state)) {
      fprintf(stderr, "ERROR. Failed getting Moon's state\n");
      return 5;
    }

    assert(std::abs(pos(0)-state(0))<1e-3);
    assert(std::abs(pos(1)-state(1))<1e-3);
    assert(std::abs(pos(2)-state(2))<1e-3);
  } /* end tests for Moon */
  
  for (int i=0; i<numtests; i++) {
    dso::MjdEpoch t = dso::MjdEpoch::random(dso::modified_julian_day(47892),
                                  dso::modified_julian_day(62502));

    if (dso::planet_pos(dso::Planet::SUN, t, pos)) {
      fprintf(stderr, "ERROR. Failed getting Sun's position\n");
      return 5;
    }

    if (dso::planet_state(dso::Planet::SUN, t, state)) {
      fprintf(stderr, "ERROR. Failed getting Sun's state\n");
      return 5;
    }

    /* Position of the Earth with respect to the Sun (GCRF). We will need this 
     * for relativistic accleration correction.
     */
    {
      double data[6];
      int pid;
      if (dso::cspice::planet_to_naif_id(dso::Planet::SUN, pid)) return 2;
      /* get position of the Earth with respect to the Sun (GCRF) [km] */
      if (dso::cspice::j2planet_pos_from(dso::cspice::mjdtt2et(t), 399, pid,
                                         data))
        return 2;
      assert(std::abs(pos(0)+data[0]*1e3)<1e-3);
      assert(std::abs(pos(1)+data[1]*1e3)<1e-3);
      assert(std::abs(pos(2)+data[2]*1e3)<1e-3);

      /* get state of Earth w.r.t Sun (GCRF) */
      double dummy;
      spkezr_c("399", dso::cspice::mjdtt2et(t), "J2000", "NONE", "10", data,
               &dummy);
      assert(std::abs(state(3)+data[3]*1e3)<1e-6);
      assert(std::abs(state(4)+data[4]*1e3)<1e-6);
      assert(std::abs(state(5)+data[5]*1e3)<1e-6);
      //printf("Note: Vearth w.r.t Sun = (%.3f, %.3f, %.3f) [m/sec]\n", data[3]*1e3, data[4]*1e3, data[5]*1e3);
      //printf("Note: Vearth w.r.t Sun = (%.3f, %.3f, %.3f) [m/sec]\n", state(3), state(4), state(5));
    }

    assert(std::abs(pos(0)-state(0))<1e-3);
    assert(std::abs(pos(1)-state(1))<1e-3);
    assert(std::abs(pos(2)-state(2))<1e-3);
  } /* end tests for Moon */

  return 0;
}
