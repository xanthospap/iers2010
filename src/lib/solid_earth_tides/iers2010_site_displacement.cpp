#include "solid_earth_tide.hpp"
#include "geodesy/core/crdtype_warppers.hpp"

Eigen::Matrix<double, 3, 1>
dso::SolidEarthTide::displacement(const dso::MjdEpoch &mjdtt,
                                  const dso::MjdEpoch &mjdut1,
                                  const Eigen::Matrix<double, 3, 1> &rsta,
                                  const Eigen::Matrix<double, 3, 1> &rMoon,
                                  const Eigen::Matrix<double, 3, 1> &rSun,
                                  const double *const delaunay_args) noexcept {
  /* Cartesian to Spherical and trig numbers */
  dso::SolidEarthTide::PointSphericalTrigs trigs_sun{
      dso::CartesianCrdConstView(rSun)};
  dso::SolidEarthTide::PointSphericalTrigs trigs_mon{
      dso::CartesianCrdConstView(rMoon)};
  dso::SolidEarthTide::PointSphericalTrigs trigs_sta{
      dso::CartesianCrdConstView(rsta)};

#ifdef DEBUG
  const auto Rt =
      dso::geodetic2lvlh(trigs_sta.msph.lat(), trigs_sta.msph.lon());
#endif
  auto dr = this->step1_displacement(rsta, rMoon, rSun, trigs_sta, trigs_mon,
                                     trigs_sun);
#ifdef DEBUG
  auto dx = Rt * dr;
  printf("STEP1         : %15.9e %15.9e %15.9e\n", dx(0), dx(1), dx(2));
#endif

  /* step-2 corrections (frequency domain) */
  dr += step2_displacement(mjdtt, mjdut1, rsta, trigs_sta, delaunay_args);
  
  /* get rotation matrix to transform between Cartesian and topocentric, i.e.
   * [x,y,z] = R *[enu]
   */
  const auto R = dso::geodetic2lvlh(trigs_sta.msph.lat(), trigs_sta.msph.lon());
  
  /* transform displacement vector from Cartesian to topocentric (enu) */
  Eigen::Matrix<double, 3, 1> dxyz = R * dr;

  /* all done */
  return dxyz;
}
