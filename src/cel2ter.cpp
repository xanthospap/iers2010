#include "cel2ter.hpp"
#include <stdexcept>

dso::Itrs2Gcrs::Itrs2Gcrs(const dso::TwoPartDate &t,
                          const dso::EopLookUpTable *eoplut)
    : eopLut(eoplut), t_tt(t) {
  if (t_tt != dso::TwoPartDate{}) {
    t_tt = dso::TwoPartDate{};
    if (prepare(t)) {
      fprintf(stderr,
              "[ERROR] Failed to construct Itrs2Gcrs instance at given date "
              "(traceback: %s)\n",
              __func__);
      throw std::runtime_error("Itrs2Gcrs() Failed ... Exiting\n");
    }
  }
}

Eigen::Matrix<double, 3, 3>
dso::detail::T(const Eigen::Matrix<double, 3, 3> &Rc2i, double era,
  const Eigen::Matrix<double, 3, 3> &Rpom) noexcept {
  // GCRS-to-ITRS matrix
  const auto Rot =
      (Rpom * (Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * Rc2i));
  return Rot;
}

/// @brief Transform a (position) vector given in ITRF to GCRF
/// @param[in] Rc2i Precession/Nutation matrix
/// @param[in] era Earth Rotation Angle (ERA) [rad]
/// @param[in] Rpom Polar motion matrix
/// @return r_gcrf = T * r_itrf
Eigen::Matrix<double, 3, 1>
dso::itrf2gcrf(const Eigen::Matrix<double, 3, 3> &Rc2i, double era,
        const Eigen::Matrix<double, 3, 3> &Rpom,
        const Eigen::Matrix<double, 3, 1> &r_itrf) noexcept {
  return detail::T(Rc2i, era, Rpom).transpose() * r_itrf;
  //const auto Rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * Rc2i;
  //return (Rpom * Rc2ti).transpose() * r_itrf;
}

Eigen::Matrix<double,3,3> dso::Itrs2Gcrs::R_gcrs2itrs() const noexcept 
{
  return detail::T(Rc2i, era, Rpom);
}

/// @brief Transform a (position) vector given in GCRF to ITRF
/// @param[in] Rc2i Precession/Nutation matrix
/// @param[in] era Earth Rotation Angle (ERA) [rad]
/// @param[in] Rpom Polar motion matrix
/// @return r_itrf = T^T * r_gcrf
Eigen::Matrix<double, 3, 1>
dso::gcrf2itrf(const Eigen::Matrix<double, 3, 3> &Rc2i, double era,
        const Eigen::Matrix<double, 3, 3> &Rpom,
        const Eigen::Matrix<double, 3, 1> &r_gcrf) noexcept {
  return detail::T(Rc2i, era, Rpom) * r_gcrf;
  // const auto Rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * Rc2i;
  // return (Rpom*Rc2ti) * r_gcrf;
}

Eigen::Matrix<double, 6, 1>
dso::gcrf2itrf(const Eigen::Matrix<double, 3, 3> &Rc2i, double era,
               const Eigen::Matrix<double, 3, 3> &Rpom, double omega_earth,
               const Eigen::Matrix<double, 6, 1> &y_gcrf) noexcept {
  // result state vector
  Eigen::Matrix<double, 6, 1> y_itrf;

  // transform position
  const Eigen::Matrix<double, 3, 1> r = y_gcrf.block<3, 1>(0, 0);
  y_itrf.block<3, 1>(0, 0) = gcrf2itrf(Rc2i, era, Rpom, r);

  // transform velocity
  const Eigen::Matrix<double, 3, 1> v = y_gcrf.block<3, 1>(3, 0);
  // instantaneous angular velocity vector
  const Eigen::Matrix<double, 3, 1> Om =
      Eigen::Matrix<double, 3, 1>({0e0, 0e0, omega_earth});
  // transformation matrix GCRS-to-TIRS
  const auto R = (Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * Rc2i);
  // transform ...
  y_itrf.block<3, 1>(3, 0) = Rpom * (R * v - Om.cross(R * r));

  return y_itrf;
}

Eigen::Matrix<double, 6, 1>
dso::itrf2gcrf(const Eigen::Matrix<double, 3, 3> &Rc2i, double era,
               const Eigen::Matrix<double, 3, 3> &Rpom, double omega_earth,
               const Eigen::Matrix<double, 6, 1> &y_itrf) noexcept {
  // result state vector
  Eigen::Matrix<double, 6, 1> y_gcrf;

  // transform position
  const Eigen::Matrix<double, 3, 1> r = y_itrf.block<3, 1>(0, 0);
  y_gcrf.block<3, 1>(0, 0) = itrf2gcrf(Rc2i, era, Rpom, r);

  // transform velocity
  const Eigen::Matrix<double, 3, 1> v = y_itrf.block<3, 1>(3, 0);
  // instantaneous angular velocity vector
  const Eigen::Matrix<double, 3, 1> Om =
      Eigen::Matrix<double, 3, 1>({0e0, 0e0, omega_earth});
  // transform ...
  y_gcrf.block<3, 1>(3, 0) =
      itrf2gcrf(Rc2i, era, Rpom, v) +
      (Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * Rc2i).transpose() *
          Om.cross(Rpom.transpose() * r);

  return y_gcrf;
}

dso::TwoPartDate dso::Itrs2Gcrs::ut1() const noexcept {
  // ∆T  = TT − UT1 = 32s.184 + ∆AT − ∆UT1 , or
  // UT1 = TAI + ΔAT - ΔUT1 , or
  // UT1 = UTC - ΔUT1
  const auto utc = t_tt.tt2tai().tai2utc().normalized();
  const auto ut1 =
      dso::TwoPartDate(utc._big, utc._small + eops.dut / dso::sec_per_day);
  return ut1.normalized();
}

int dso::Itrs2Gcrs::prepare(const dso::TwoPartDate &tt_mjd) noexcept {
  if (tt_mjd != t_tt) {
    // interpolate EOPs, using a Lagrange polynomial of 5th degree and
    // ASSUMING that the EOP table instance is tabulated in TT
    if (eopLut->interpolate(tt_mjd, eops, 5, true)) {
      fprintf(
          stderr,
          "[ERROR] Failed interpolating EOPs for given date (traceback: %s)\n",
          __func__);
      return 1;
    }
    //printf("Interpolated eop: %.4f %.4f %.6f %.6f %.4f %.4f %.5e\n", eops.xp, eops.yp, eops.dut, eops.lod, eops.dx, eops.dy, eops.omega());

    // set ERA(t) angle
    {
      // temporarily set instance's date to compute ERA, them swap back
      const auto t_temp = t_tt;
      t_tt = tt_mjd;
      era = iers2010::sofa::era00(this->ut1());
      t_tt = t_temp;
    }

    // Construct Q(t) matrix (precession/nutation)
    {
      // call the routine XY06 to obtain the IAU 2006/2000A X, Y from series,
      // i.e. CIP coordinates
      double X, Y;
      iers2010::sofa::xy06(tt_mjd, X, Y);

      // call the routine S06 to obtain s
      const double s = iers2010::sofa::s06(tt_mjd, X, Y);

      // apply CIP corrections, ΔX and ΔY (if any)
      X += dso::sec2rad(eops.dx);
      Y += dso::sec2rad(eops.dy);

      // call routine C2IXYS, giving the GCRS-to-CIRS matrix Q
      Rc2i = iers2010::sofa::c2ixys(X, Y, s);
    }

    // construct polar motion matrix W(t)
    {
      const double sp = iers2010::sofa::sp00(tt_mjd);
      Rpom = iers2010::sofa::pom00(dso::sec2rad(eops.xp), dso::sec2rad(eops.yp),
                                sp);
    }

    // set current time
    t_tt = tt_mjd;
  }
  return 0;
}
