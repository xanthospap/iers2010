#include "iau.hpp"

void iers2010::sofa::pn00(const dso::TwoPartDate &mjd_tt, double dpsi,
                          double deps, double &epsa,
                          Eigen::Matrix<double, 3, 3> &rb,
                          Eigen::Matrix<double, 3, 3> &rp,
                          Eigen::Matrix<double, 3, 3> &rbp,
                          Eigen::Matrix<double, 3, 3> &rn,
                          Eigen::Matrix<double, 3, 3> &rbpn) noexcept {
  // IAU 2000 precession-rate adjustments.
  double dpsipr, depspr;
  iers2010::sofa::pr00(mjd_tt, dpsipr, depspr);

  // Mean obliquity, consistent with IAU 2000 precession-nutation.
  epsa = iers2010::sofa::obl80(mjd_tt) + depspr;

  // Frame bias and precession matrices and their product.
  iers2010::sofa::bp00(mjd_tt, rb, rp, rbp);

  // Nutation matrix.
  rn = iers2010::sofa::numat(epsa, dpsi, deps);

  // Bias-precession-nutation matrix (classical).
  rbpn = rn * rbp;
}
