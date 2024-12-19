#include "planets.hpp"

int dso::planet_pos(dso::Planet p, const dso::MjdEpoch &mjd_tt,
                    Eigen::Matrix<double, 3, 1> &pos) noexcept {
  double data[3];
  int pid;

  if (cspice::planet_to_naif_id(p, pid))
    return 2;

  /* get position in [km] */
  int status =
      cspice::j2planet_pos_from(cspice::mjdtt2et(mjd_tt), pid, 399, data);

  /* assign to output matrix */
  pos(0) = data[0] * 1e3;
  pos(1) = data[1] * 1e3;
  pos(2) = data[2] * 1e3;

  /* return */
  return status;
}
