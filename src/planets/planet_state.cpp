#include "planets.hpp"

int dso::planet_state(dso::Planet p, const dso::MjdEpoch &mjd_tt,
                    Eigen::Matrix<double, 6, 1> &pv) noexcept {
  double data[6];
  
  /* need to transform the target/observer to chars, holding their (int) id's */
  char trg[8];
  if (int errc = dso::cspice::planet_to_naif_idstr(p, trg, 8); errc) {
    fprintf(stderr,
            "[ERROR] Failed matching NAIF id for given planet; error=%d! "
            "(traceback: %s)\n",
            errc, __func__);
    return 1;
  }

  /* get state in [km] and [km/sec] */
  int status =
      cspice::j2planet_state_from(cspice::mjdtt2et(mjd_tt), trg, "399", data);

  /* assign to output matrix */
  pv(0) = data[0] * 1e3;
  pv(1) = data[1] * 1e3;
  pv(2) = data[2] * 1e3;
  pv(3) = data[3] * 1e3;
  pv(4) = data[4] * 1e3;
  pv(5) = data[5] * 1e3;

  /* return */
  return status;
}
