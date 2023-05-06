#include "iau.hpp"

void iers2010::sofa::bp00(const dso::TwoPartDate &mjd_tt,
                          Eigen::Matrix<double, 3, 3> &rb,
                          Eigen::Matrix<double, 3, 3> &rp,
                          Eigen::Matrix<double, 3, 3> &rbp) noexcept {
  /* J2000.0 obliquity (Lieske et al. 1977) */
  constexpr double EPS0 = 84'381.448e0 * iers2010::DAS2R;

  /* Interval between fundamental epoch J2000.0 and current date (JC) */
  const double t = mjd_tt.jcenturies_sinceJ2000();

  /* Frame bias */
  double dpsibi, depsbi, dra0;
  iers2010::sofa::bi00(dpsibi, depsbi, dra0);

  /* Precession angles (Lieske et al. 1977) */
  const double psia77 = (5'038.7784e0 + (-1.07259e0 + (-0.001147e0) * t) * t) *
                        t * iers2010::DAS2R;
  const double oma77 =
      EPS0 + ((0.05127e0 + (-0.007726e0) * t) * t) * t * iers2010::DAS2R;
  const double chia =
      (10.5526e0 + (-2.38064e0 + (-0.001125e0) * t) * t) * t * iers2010::DAS2R;

  /* Apply IAU 2000 precession corrections */
  double dpsipr, depspr;
  iers2010::sofa::pr00(mjd_tt, dpsipr, depspr);
  const double psia = psia77 + dpsipr;
  const double oma = oma77 + depspr;

  /* Frame bias matrix: GCRS to J2000.0 */
  rb = Eigen::AngleAxisd(depsbi, Eigen::Vector3d::UnitX()) *
       Eigen::AngleAxisd(-(dpsibi * std::sin(EPS0)), Eigen::Vector3d::UnitY()) *
       Eigen::AngleAxisd(-dra0, Eigen::Vector3d::UnitZ());

  /* Precession matrix: J2000.0 to mean of date. */
  rp = Eigen::AngleAxisd(-chia, Eigen::Vector3d::UnitZ()) *
       Eigen::AngleAxisd(oma, Eigen::Vector3d::UnitX()) *
       Eigen::AngleAxisd(psia, Eigen::Vector3d::UnitZ()) *
       Eigen::AngleAxisd(-EPS0, Eigen::Vector3d::UnitX());

  /* Bias-precession matrix: GCRS to mean of date. */
  rbp = rp * rb;
}
