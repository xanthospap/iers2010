#include "solid_earth_tide.hpp"
#include <array>

int dso::SolidEarthTide::stokes_coeffs(
    const dso::MjdEpoch &mjdtt, const dso::MjdEpoch &mjdut1,
    const Eigen::Matrix<double, 3, 1> &rMoon,
    const Eigen::Matrix<double, 3, 1> &rSun,
    const double *const delaunay_args) noexcept {
  /* store dC and dS here, from Step-1 */
  std::array<double, 12> dC, dS;

  /* compute Step-1 corrections (Sun+Moon) */
  potential_step1(rMoon, rSun, dC, dS);

  /* compute Step-2 corrections */
  double dC20, dC21, dS21, dC22, dS22;
  potential_step2(mjdtt, mjdut1, delaunay_args, dC20, dC21, dS21, dC22,
                         dS22);

  /* assign corrections to the isntance's Stokes coefficients */
  mcs.clear();
  /* C coeffs */
  mcs.C(2, 0) = dC[0] + dC20;
  mcs.C(3, 0) = dC[3];
  mcs.C(4, 0) = dC[7];
  mcs.C(2, 1) = dC[1] + dC21;
  mcs.C(3, 1) = dC[4];
  mcs.C(4, 1) = dC[8];
  mcs.C(2, 2) = dC[2] + dC22;
  mcs.C(3, 2) = dC[5];
  mcs.C(4, 2) = dC[9];
  mcs.C(3, 3) = dC[6];
  mcs.C(4, 3) = dC[10];
  mcs.C(4, 4) = dC[11];

  /* S coeffs */
  mcs.S(2, 1) = dS[1] + dS21;
  mcs.S(3, 1) = dS[4];
  mcs.S(4, 1) = dS[8];
  mcs.S(2, 2) = dS[2] + dS22;
  mcs.S(3, 2) = dS[5];
  mcs.S(4, 2) = dS[9];
  mcs.S(3, 3) = dS[6];
  mcs.S(4, 3) = dS[10];
  mcs.S(4, 4) = dS[11];


  //printf("C(2,0)=%+.15f C(2,1)=%+.15f C(2,2)=%+.15f\n", mcs.C(2,0), mcs.C(2,1), mcs.C(2,2));
  //printf("C(3,0)=%+.15f C(3,1)=%+.15f C(3,2)=%+.15f C(3,3)=%+.15f\n", mcs.C(3,0), mcs.C(3,1), mcs.C(3,2), mcs.C(3,3));
  //printf("C(4,0)=%+.15f C(4,1)=%+.15f C(4,2)=%+.15f C(4,3)=%+.15f C(4,4)=%+.15f\n", mcs.C(4,0), mcs.C(4,1), mcs.C(4,2), mcs.C(4,3), mcs.C(4,4));

  return 0;
}
