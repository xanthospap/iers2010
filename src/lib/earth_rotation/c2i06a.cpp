#include "earth_rotation.hpp"
#include "iau.hpp"

Eigen::Quaterniond dso::c2i06a(const dso::MjdEpoch &tt,
                               const dso::EopRecord &eops) noexcept {
  double fargs[14];
  double Xcip, Ycip;

  /* compute (X,Y) CIP and fundamental arguments */
  dso::xycip06a(tt, Xcip, Ycip, fargs);

  /* use fundamental arguments to compute s */
  const double s = dso::s06(tt, Xcip, Ycip, fargs);

  /* apply CIP corrections */
  Xcip += dso::sec2rad(eops.dX());
  Ycip += dso::sec2rad(eops.dY());

  return dso::detail::W(dso::sec2rad(eops.xp()), dso::sec2rad(eops.yp()),
                        dso::sp00(tt)) *
         dso::detail::R(dso::era00(tt.tt2ut1(eops.dut()))) *
         dso::detail::C(Xcip, Ycip, s);
}

Eigen::Quaterniond dso::detail::c2i(double era, double s, double sp,
                                    double Xcip, double Ycip, double xp,
                                    double yp) noexcept {
  /* compute full rotation matrix: [C x R x W] */
  return dso::detail::W(xp, yp, sp) * dso::detail::R(era) *
         dso::detail::C(Xcip, Ycip, s);
}

Eigen::Quaterniond
dso::detail::c2i(double era, double s, double sp, double Xcip, double Ycip,
                 double xp, double yp, double lod,
                 Eigen::Matrix<double, 3, 3> &dRdt) noexcept {
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

    /*
   Eigen::Matrix<double, 3, 3> Rdot = dso::detail::R(era).toRotationMatrix() *
                                      (dso::earth_rotation_rate(lod));
   const double pcos = Rdot(0, 0);
   const double msin = Rdot(0, 1);
   Rdot(0, 0) = -msin;
   Rdot(1, 1) = -msin;
   Rdot(0, 1) = -pcos;
   Rdot(1, 0) = pcos;
   Rdot(2, 2) = 0e0;
   */
    Eigen::Matrix<double, 3, 3> S = Eigen::Matrix<double, 3, 3>::Zero();
    S(0, 1) = 1e0;
    S(1, 0) = -1e0;
    const Eigen::Matrix<double,3,3> M = dso::earth_rotation_rate(lod) * S *
                   dso::detail::R(era).toRotationMatrix();
    printf("M=\n");
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        printf("+%.3e ", M(i, j));
      }
      printf("\n");
    }

    dRdt = dso::detail::W(xp, yp, sp).toRotationMatrix() * M *
           dso::detail::C(Xcip, Ycip, s).toRotationMatrix();
  }

  return dso::detail::c2i(era, s, sp, Xcip, Ycip, xp, yp);
}

Eigen::Quaterniond dso::c2i06a(const dso::MjdEpoch &tt,
                               const dso::EopRecord &eops,
                               Eigen::Matrix<double, 3, 3> &dRdt) noexcept {
  double fargs[14];
  double Xcip, Ycip;

  /* compute (X,Y) CIP and fundamental arguments */
  dso::xycip06a(tt, Xcip, Ycip, fargs);

  /* use fundamental arguments to compute s */
  const double s = dso::s06(tt, Xcip, Ycip, fargs);

  /* apply CIP corrections */
  Xcip += dso::sec2rad(eops.dX());
  Ycip += dso::sec2rad(eops.dY());

  return dso::detail::c2i(dso::era00(tt.tt2ut1(eops.dut())), s, dso::sp00(tt),
                          Xcip, Ycip, dso::sec2rad(eops.xp()),
                          dso::sec2rad(eops.yp()), eops.lod(), dRdt);
}
