#include "eop.hpp"
#include "iau.hpp"
#include "datetime/dtcalendar.hpp"
#include <stdexcept>

/*
 * According to IERS 2010, the transformation to be used to relate the 
 * International Terrestrial Reference System (ITRS) to the Geocentric 
 * Celestial Reference System (GCRS) at the date t of the observation can be 
 * written as:
 *              [GCRS] = Q(t) R(t) W(t) [ITRS]                           (5.1)
 * where Q(t), R(t) and W (t) are the transformation matrices arising from the 
 * motion of the celestial pole in the celestial reference system, from the 
 * rotation of the Earth around the axis associated with the pole, and from 
 * polar motion, respectively.
 * 
 * The parameter t, used in Eq. (5.1) as well as in the following expressions, 
 * is defined by:
 *              t = (TT - 2000 January 1d 12h TT) in days / 36525        (5.2)
 *
 * Matrix W (polar motion matrix)
 * ------------------------------
 * According to IAU 2006 Resolution B2, the system at date t as realized from 
 * the ITRS by applying the transformation W(t) in both procedures is the 
 * “Terrestrial Intermediate Reference System” (TIRS). It uses the CIP as its 
 * z-axis and the TIO as its x-axis.
 * The transformation matrix arising from polar motion (i.e. relating ITRS and 
 * TIRS) can be expressed as:
 *              W(t) = R3(−s0 ) · R2(xp) · R1(yp)                        (5.3)
 * Note that the function rpom00 actually computes the matrix W^T, so that
 *              [ITRS] = RPOM [TIRS]
 * 
 * Matrix R (Earth Rotation)
 * -------------------------
 * The CIO based transformation matrix arising from the rotation of the Earth 
 * around the axis of the CIP (i.e. relating TIRS and CIRS), can be expressed 
 * as:
 *              R(t) = R3 (−ERA)                                         (5.5)
 * where ERA is the Earth Rotation Angle between the CIO and the TIO at date t 
 * on the equator of the CIP, which is the rigorous definition of the sidereal 
 * rotation of the Earth.
 *
 * Matrix Q (precession/nutation matrix)
 * -------------------------------------
 * The CIO based procedure realizes an intermediate celestial reference system 
 * at date t that uses the CIP as its z-axis and the CIO as its x-axis. 
 * According to IAU 2006 Resolution B2, it is called the “Celestial 
 * Intermediate Reference System” (CIRS). It uses the Earth Rotation Angle in 
 * the transformation matrix R(t), and the two coordinates of the CIP in the 
 * GCRS in the transformation matrix Q(t). 
 * The CIO based transformation matrix arising from the motion of the CIP in 
 * the GCRS (i.e. relating CIRS and GCRS), can be expressed as:
 *              Q(t) = R3(−E) · R2(−d) · R3(E) · R3(s)                   (5.6)
 * Note that the function c2ixys the matrix Q^T, so that
 *              [ITRS] = RPOM * R3(ERA) RC2I [GCRS]
 * 
 * 
 * Notes:
 *  All (passed around) datetime instances in here are considered to be in 
 *  the TT time-scale, except if explicitelly otherwise stated.
 * 
 *  Each of the transformation matrix components W(t), R(t) and Q(t) of 
 *  Eq.(5.1) is a series of rotations about the axes 1, 2 and 3 of the 
 *  coordinate frame. In the following, R1 , R2 and R3 denote rotation 
 *  matrices with positive angle about the axes 1, 2 and 3. 
 */

namespace dso {

namespace detail {
/// @brief Tranformation matrix to be used for ITRF to GCRF transformation.
/// So that (according to IERS-2010), the relationship for transforming a
/// position vector between the two frames, is:
///      r_gcrf = Q(t) * R(t) * W(t) * r_itrf,
/// where 
///  * Q is the matrix to model the motion of the celestial pole in the 
///    celestial reference system
///  * R transformation matrix from the rotation of the Earth around the axis 
///    associated with the pole, that R <- R_z(-ERA)
///  * W transformation matrix to account for the polar motion
/// This function will compute the matrix:
///      T <- Q(t) * R(t) * W(t)
///
/// @param[in] Q Precession/Nutation matrix
/// @param[in] era Earth Rotation Angle (ERA) [rad]
/// @param[in] W Polar motion matrix
/// @return T = Q * R * W
Eigen::Matrix<double, 3, 3>
T(const Eigen::Matrix<double, 3, 3> &Q, double era,
  const Eigen::Matrix<double, 3, 3> &W) noexcept;
} // namespace detail

/// @brief Transform a (position) vector given in ITRF to GCRF
/// @param[in] Q Precession/Nutation matrix
/// @param[in] era Earth Rotation Angle (ERA) [rad]
/// @param[in] W Polar motion matrix
/// @return r_gcrf = T * r_itrf
Eigen::Matrix<double, 3, 1>
itrf2gcrf(const Eigen::Matrix<double, 3, 3> &Q, double era,
        const Eigen::Matrix<double, 3, 3> &W,
        const Eigen::Matrix<double, 3, 1> &r_itrf) noexcept;

/// @brief Transform a (position) vector given in GCRF to ITRF
/// @param[in] Q Precession/Nutation matrix
/// @param[in] era Earth Rotation Angle (ERA) [rad]
/// @param[in] W Polar motion matrix
/// @return r_itrf = T^T * r_gcrf
Eigen::Matrix<double, 3, 1>
gcrf2itrf(const Eigen::Matrix<double, 3, 3> &Q, double era,
        const Eigen::Matrix<double, 3, 3> &W,
        const Eigen::Matrix<double, 3, 1> &r_gcrf) noexcept;

/// @brief Transform a state (position/velocity) vector given in GCRF to ITRF
/// @param[in] Q Precession/Nutation matrix
/// @param[in] era Earth Rotation Angle (ERA) [rad]
/// @param[in] W Polar motion matrix
/// @param[in] omega_earth Instantaneous rotation rate of Earth [rad/sec]
/// @return A 6-element vector [r_itrf, v_itrf], where 
///         r_itrf = T^T * r_gcrf
///         v_itrf = W^T [R^T Q^T v_gcrf - ω x r_tirs]
/// See Vallado 2013, Section 3.7, "Implementing the IAU-2010 Conventions"
Eigen::Matrix<double, 6, 1>
gcrf2itrf(const Eigen::Matrix<double, 3, 3> &Q, double era,
        const Eigen::Matrix<double, 3, 3> &W, double omega_earth,
        const Eigen::Matrix<double, 6, 1> &y_gcrf) noexcept;

/// @brief Transform a state (position/velocity) vector given in ITRF to GCRF
/// @param[in] Q Precession/Nutation matrix
/// @param[in] era Earth Rotation Angle (ERA) [rad]
/// @param[in] W Polar motion matrix
/// @param[in] omega_earth Instantaneous rotation rate of Earth [rad/sec]
/// @return A 6-element vector [r_gcrf, v_gcrf], where 
///         r_gcrf = T * r_itrf
///         v_gcrf = Q * R [W v_itrf + ω x r_tirs]
/// See Vallado 2013, Section 3.7, "Implementing the IAU-2010 Conventions"
Eigen::Matrix<double, 6, 1>
itrf2gcrf(const Eigen::Matrix<double, 3, 3> &Q, double era,
        const Eigen::Matrix<double, 3, 3> &W, double omega_earth,
        const Eigen::Matrix<double, 6, 1> &y_itrf) noexcept;

class Itrs2Gcrs {
  const dso::EopLookUpTable * const eopLut;
  ///< Current epoch in TT (MJD)
  dso::TwoPartDate t_tt;
  ///< Current EOPs/ERPs after interpolation
  dso::EopRecord eops;
  ///< Polar motion matrix: transforms TIRS-to-ITRS (i.e. W^T in 5.1)
  Eigen::Matrix<double,3,3> Rpom;
  ///< Precession/Nutation matrix: transforms GCRS-to-CIRS matrix (i.e. Q^t in 5.1)
  Eigen::Matrix<double,3,3> Rc2i;
  ///< Eart rotation angle in [rad]. (R3(era) transforms CIRS-to-TIRS)
  double era;

public:
  Itrs2Gcrs(const dso::TwoPartDate &t = dso::TwoPartDate{},
            const dso::EopLookUpTable *eoplut = nullptr);

  dso::TwoPartDate ut1() const noexcept;
  int prepare(const dso::TwoPartDate &tt_mjd) noexcept;
  dso::EopRecord eop() const noexcept {return eops;}
  Eigen::Matrix<double,3,3> gcrf2itrf() const noexcept;
  Eigen::Matrix<double,3,3> gcrf2tirs() const noexcept;
  Eigen::Matrix<double,3,3> itrf2gcrf() const noexcept;
  const Eigen::Matrix<double,3,3> &rpom() const noexcept {return Rpom;}
  const Eigen::Matrix<double,3,3> &rc2i() const noexcept {return Rc2i;}
  double earth_rotation_angle() const noexcept {return era;}
  /// @return Angular velocity of Earth (Ω) in [rad/sec]
  double omega_earth() const noexcept {return eops.omega();}

  Eigen::Matrix<double,3,3> ddt_gcrf2itrf() const noexcept;

  Eigen::Matrix<double, 3, 1>
  itrf2gcrf(const Eigen::Matrix<double, 3, 1> &ritrf) const noexcept;
  Eigen::Matrix<double, 6, 1>
  itrf2gcrf(const Eigen::Matrix<double, 6, 1> &yitrf) const noexcept;
  Eigen::Matrix<double, 3, 1>
  gcrf2itrf(const Eigen::Matrix<double, 3, 1> &rgcrf) const noexcept;
  Eigen::Matrix<double, 6, 1>
  gcrf2itrf(const Eigen::Matrix<double, 6, 1> &ygcrf) const noexcept;
}; // class Itrs2Gcrs

} // namespace dso
