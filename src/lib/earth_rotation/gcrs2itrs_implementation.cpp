#include "earth_rotation.hpp"
#include "iau.hpp"

Eigen::Matrix<double, 3, 3> dso::detail::gcrs2itrs(double era, double s,
                                                   double sp, double Xcip,
                                                   double Ycip, double xp,
                                                   double yp) noexcept {
  /* compute full rotation matrix: [C x R x W] */
  auto Rm = W(xp, yp, sp) * detail::R(era) * detail::C(Xcip, Ycip, s);

  return Rm.toRotationMatrix();
}

Eigen::Matrix<double, 3, 3>
dso::detail::gcrs2itrs(double era, double s, double sp, double Xcip,
                       double Ycip, double xp, double yp, double lod,
                       Eigen::Matrix<double, 3, 3> &dRdt) noexcept {
  /* compute full rotation matrix: [C x R x W] */
  auto Rm = detail::W(xp, yp, sp) * detail::R(era) * detail::C(Xcip, Ycip, s);

  // Eigen::Matrix<double,3,3> T = Eigen::Matrix<double,3,3>::Zero();
  // T(0,1) = 1e0;
  // T(1,0) = -1e0;
  // dRdt = W(xp, yp, sp) * (earth_rotation_rate(lod)*T*R(era)) * C(Xcip, Ycip,
  // s);

  {
    /* compute simplified rotation matrix derivative, where we only consider
     * dR/dt: i.e.
     * [C x dR3/dt x W], where dR3/dt =
     *
     * | -sin(era) -cos(era) 0 |
     * |  cos(era) -sin(era) 0 | * ω_{earth}
     * |    0         0      0 |
     *
     * and ω_{earth} = earth_rotation_rate(lod) in [rad/sec]
     */
    Eigen::Matrix<double, 3, 3> Rdot =
        detail::R(era).toRotationMatrix() * (-earth_rotation_rate(lod));
    const double pcos = Rdot(0, 0);
    const double msin = Rdot(0, 1);
    Rdot(0, 0) = msin;
    Rdot(1, 1) = msin;
    Rdot(0, 1) = -pcos;
    Rdot(1, 0) = pcos;
    Rdot(2, 2) = 0e0;
    dRdt = detail::W(xp, yp, sp) * Rdot * detail::C(Xcip, Ycip, s);
  }

  return Rm.toRotationMatrix();
}

Eigen::Matrix<double, 3, 3>
dso::detail::gcrs2tirs(double era, double s, double sp, double Xcip,
                       double Ycip, double xp, double yp,
                       Eigen::Matrix<double, 3, 3> *Wpom) noexcept {
  /* compute full rotation matrix: [C x R] */
  auto Rm = detail::R(era) * detail::C(Xcip, Ycip, s);

  /* polar motion matrix */
  if (Wpom)
    *Wpom = detail::W(xp, yp, sp);

  return Rm.toRotationMatrix();
}
