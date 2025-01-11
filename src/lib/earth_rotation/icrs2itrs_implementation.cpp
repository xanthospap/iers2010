#include "earth_rotation.hpp"
#include "iau.hpp"

namespace {

/* polar motion matrix */
auto W(double xp, double yp, double sp) noexcept {
  /* W = R3(-s') R2(xp) R1(yp) */
  return (Eigen::AngleAxisd(yp, Eigen::Vector3d::UnitX()) * 
    Eigen::AngleAxisd(xp, Eigen::Vector3d::UnitY()) * 
    Eigen::AngleAxisd(-sp, Eigen::Vector3d::UnitZ()));
}

auto C(double x, double y, double s) noexcept {
  /* Obtain the spherical angles E and d:
   * X = sin d cos E,
   * Y = sin d sin E,
   * Z = cos d
   */
  const double r2 = x * x + y * y;
  const double e = (r2 > 0e0) ? std::atan2(y, x) : 0e0;
  const double d = std::atan(std::sqrt(r2 / (1e0 - r2)));

  /* R3 (−E) x R2 (−d) x R3 (E) x R3 (s) */
  return (Eigen::AngleAxisd(e + s, Eigen::Vector3d::UnitZ()) *
          Eigen::AngleAxisd(-d, Eigen::Vector3d::UnitY()) *
          Eigen::AngleAxisd(-e, Eigen::Vector3d::UnitZ()));
}

auto R(double era) noexcept {
  return Eigen::AngleAxisd(-era, Eigen::Vector3d::UnitZ());
}

} /* unnamed namespace */

Eigen::Matrix<double, 3, 3> dso::detail::gcrs2itrs(double era, double s,
                                                   double sp, double Xcip,
                                                   double Ycip, double xp,
                                                   double yp) noexcept {
  /* compute full rotation matrix: [C x R x W] */
  auto Rm =  W(xp, yp, sp) * R(era) * C(Xcip, Ycip, s);

  return Rm.toRotationMatrix();
}

Eigen::Matrix<double, 3, 3>
dso::detail::gcrs2itrs(double era, double s, double sp, double Xcip,
                       double Ycip, double xp, double yp, double lod,
                       Eigen::Matrix<double, 3, 3> &dRdt) noexcept {
  /* compute full rotation matrix: [C x R x W] */
  auto Rm = W(xp, yp, sp) * R(era) * C(Xcip, Ycip, s);

  Eigen::Matrix<double,3,3> T = Eigen::Matrix<double,3,3>::Zero();
  T(0,1) = 1e0;
  T(1,0) = -1e0;
  dRdt = W(xp, yp, sp) * (earth_rotation_rate(lod)*T*R(era)) * C(Xcip, Ycip, s);

  //{
  //  /* compute simplified rotation matrix derivative, where we only consider
  //   * dR/dt: i.e.
  //   * [C x dR3/dt x W], where dR3/dt =
  //   *
  //   * | -sin(era) -cos(era) 0 |
  //   * |  cos(era) -sin(era) 0 | * ω_{earth}
  //   * |    0         0      0 |
  //   *
  //   * and ω_{earth} = earth_rotation_rate(lod) in [rad/sec]
  //   */
  //  Eigen::Matrix<double, 3, 3> Rdot =
  //      R(era).toRotationMatrix() * earth_rotation_rate(lod);
  //  const double pcos = Rdot(0, 0);
  //  const double msin = Rdot(0, 1);
  //  Rdot(0, 0) = msin;
  //  Rdot(1, 1) = msin;
  //  Rdot(0, 1) = -pcos;
  //  Rdot(1, 0) = pcos;
  //  Rdot(2, 2) = 0e0;
  //  dRdt = W(xp, yp, sp) * Rdot * C(Xcip, Ycip, s);
  //}

  return Rm.toRotationMatrix();
}

Eigen::Matrix<double, 3, 3>
dso::detail::gcrs2tirs(double era, double s, double sp, double Xcip,
                       double Ycip, double xp, double yp,
                       Eigen::Matrix<double, 3, 3> *Wpom) noexcept {
  /* compute full rotation matrix: [C x R] */
  auto Rm = R(era) * C(Xcip, Ycip, s);

  /* polar motion matrix */
  if (Wpom) *Wpom = W(xp, yp, sp);

  return Rm.toRotationMatrix();
}
