#include "iau.hpp"

void iers2010::sofa::pn06(const dso::TwoPartDate &mjd_tt, double dpsi,
                          double deps, double &epsa,
                          Eigen::Matrix<double, 3, 3> &rb,
                          Eigen::Matrix<double, 3, 3> &rp,
                          Eigen::Matrix<double, 3, 3> &rbp,
                          Eigen::Matrix<double, 3, 3> &rn,
                          Eigen::Matrix<double, 3, 3> &rbpn) noexcept {
  // Bias-precession Fukushima-Williams angles of J2000.0 = frame bias.
  double gamb, phib, psib, eps;
  //iers2010::sofa::pfw06(iers2010::DJM0, iers2010::DJM00, gamb, phib, psib, eps);
  iers2010::sofa::pfw06(dso::TwoPartDate(dso::j2000_mjd,0e0), gamb, phib, psib, eps);

  // B matrix.
  rb = iers2010::sofa::fw2m(gamb, phib, psib, eps);

  // Bias-precession Fukushima-Williams angles of date.
  iers2010::sofa::pfw06(mjd_tt, gamb, phib, psib, eps);

  // Bias-precession matrix.
  rbp = iers2010::sofa::fw2m(gamb, phib, psib, eps);

  // Solve for precession matrix.
  rp = rbp * rb.transpose();

  // Equinox-based bias-precession-nutation matrix.
  rbpn = iers2010::sofa::fw2m(gamb, phib, psib + dpsi, eps + deps);

  // Solve for nutation matrix.
  rn = rbpn * rbp.transpose();

  // Obliquity, mean of date.
  epsa = eps;

  // Finished.
  return;
}
