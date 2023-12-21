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
  solid_earth_tide_step1(rMoon, rSun, dC, dS);

  /* compute Step-2 corrections */
  double dC20, dC21, dS21, dC22, dS22;
  solid_earth_tide_step2(mjdtt, mjdut1, delaunay_args, dC20, dC21, dS21, dC22,
                         dS22);

  /* assign corrections to the isntance's Stokes coefficients */
  m_cs.clear();
  /* C coeefs */
  m_cs.C(2, 0) = dC[0] + dC20;
  m_cs.C(3, 0) = dC[3];
  m_cs.C(4, 0) = dC[7];
  m_cs.C(2, 1) = dC[1] + dC21;
  m_cs.C(3, 1) = dC[4];
  m_cs.C(4, 1) = dC[8];
  m_cs.C(2, 2) = dC[2] + dC22;
  m_cs.C(3, 2) = dC[5];
  m_cs.C(4, 2) = dC[9];
  m_cs.C(3, 3) = dC[6];
  m_cs.C(4, 3) = dC[10];
  m_cs.C(4, 4) = dC[11];

  /* S coeefs */
  m_cs.S(2, 1) = dS[1] + dS21;
  m_cs.S(3, 1) = dS[4];
  m_cs.S(4, 1) = dS[8];
  m_cs.S(2, 2) = dS[2] + dS22;
  m_cs.S(3, 2) = dS[5];
  m_cs.S(4, 2) = dS[9];
  m_cs.S(3, 3) = dS[6];
  m_cs.S(4, 3) = dS[10];
  m_cs.S(4, 4) = dS[11];

  return 0;
}
