#include "iau.hpp"

Eigen::Quaterniond dso::detail::itrs2gcrs_quaternion(double era, double s,
                                                     double sp, double Xcip,
                                                     double Ycip, double xp,
                                                     double yp) noexcept {

  /* Equation (14) */
  const double cxp2 = std::cos(xp / 2e0);
  const double sxp2 = std::sin(xp / 2e0);
  const double cyp2 = std::cos(yp / 2e0);
  const double syp2 = std::sin(yp / 2e0);
  const double tw = cxp2 * cyp2;
  const double aw = cxp2 * syp2;
  const double bw = cyp2 * sxp2;
  const double cw = sxp2 * syp2;

  /* compute θ' = -(s - θ -s') = θ + s' - s, Equation (10) */
  const double thetap = era + sp - s;
  const double ct2 = std::cos(thetap / 2e0);
  const double st2 = std::sin(thetap / 2e0);

  /* Obtain the spherical angles E and d, Equation (5) */
  const double r2 = Xcip * Xcip + Ycip * Ycip;
  const double e = (r2 > 0e0) ? std::atan2(Ycip, Xcip) : 0e0;
  const double d = std::atan(std::sqrt(r2 / (1e0 - r2)));
  const double X = std::sin(d) * std::cos(e);
  const double Y = std::sin(d) * std::sin(e);
  const double Z = std::cos(d);

  /* Components of quaternion, Equation (13) */
  const double _1pZ = (1e0 + Z);
  const double fac = 1e0 / (std::sqrt(2e0 * _1pZ));
  const double T = fac * (ct2 * (X * bw - Y * aw + _1pZ * tw) +
                          st2 * (X * aw + Y * bw + _1pZ * cw));
  const double A = fac * (ct2 * (X * cw + Y * tw + _1pZ * aw) +
                          st2 * (-X * tw + Y * cw - _1pZ * bw));
  const double B = fac * (ct2 * (-X * tw + Y * cw + _1pZ * bw) +
                          st2 * (-X * cw - Y * tw + _1pZ * aw));
  const double C = fac * (ct2 * (-X * aw - Y * bw + _1pZ * cw) +
                          st2 * (X * bw - Y * aw - _1pZ * tw));

  /* return quaternion */
  return Eigen::Quaterniond(T, A, B, C);
}
