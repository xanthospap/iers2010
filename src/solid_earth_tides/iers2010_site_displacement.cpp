#include "solid_earth_tide.hpp"

Eigen::Matrix<double, 3, 1>
dso::SolidEarthTide::displacement(const dso::MjdEpoch &mjdtt,
                                  const dso::MjdEpoch &mjdut1,
                                  const Eigen::Matrix<double, 3, 1> &rsta,
                                  const Eigen::Matrix<double, 3, 1> &rMoon,
                                  const Eigen::Matrix<double, 3, 1> &rSun,
                                  const double *const delaunay_args) noexcept {
  /* Cartesian to Spherical and trig numbers */
  dso::SolidEarthTide::PointSphericalTrigs trigs_sun(rSun);
  dso::SolidEarthTide::PointSphericalTrigs trigs_mon(rMoon);
  dso::SolidEarthTide::PointSphericalTrigs trigs_sta(rsta);

  auto dr = this->step1_displacement(rsta, rMoon, rSun, trigs_sta, trigs_mon,
                                     trigs_sun);

  /* step-2 corrections (frequency domain) */
  dr += step2_displacement(mjdtt, mjdut1, rsta, trigs_sta, delaunay_args);

  /* all done */
  return dr;
}
