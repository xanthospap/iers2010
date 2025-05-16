/** @file
 *
 * Class and functions to implement an easy-to-use API to describe Earth
 * rotation and the ties between ITRS and GCRS frames.
 */

#ifndef __DSO_EARTH_ROTATION_API_HPP__
#define __DSO_EARTH_ROTATION_API_HPP__

#include "eigen3/Eigen/Eigen"
#include "eigen3/Eigen/src/Geometry/Quaternion.h"
#include "eop.hpp"
#include "iersconst.hpp"

namespace dso {

/** @brief Celestial-to-terrestrial transformation, i.e. [TRS] = q * [CRS]
 *
 * The transformation follows the IAU 2006/2000A precession-nutation model.
 *
 * @param[in] tt Epoch of request, in TT scale.
 * @param[in] eop An EopRecord instance, which will be used to extract: pole
 *            motion, i.e. xp, yp, corrections to CIP, i.e. dX and dY and
 *            Delta UT1.
 * @return The rotation quaternion to transform a vector in the sense:
 *         [TRS] = q * [CRS]
 */
Eigen::Quaterniond c2i06a(const MjdEpoch &tt, const EopRecord &eop) noexcept;

/** @brief Celestial-to-terrestrial transformation for state, 
 * i.e. v[TRS] = q * v[CRS] + dM * r[CRS]
 *
 * The transformation follows the IAU 2006/2000A precession-nutation model.
 * The matrix dM (i.e. dMdt) is the derivative of M, where M operates in the 
 * sense: r[TRS] = M * r[CRS]. For the derivative, we only consider Earth 
 * rotation, with d(ERA)/dt = omega. Hence, if M is Q*R*W (sticking to IERS 
 * 2010 notation), with R = R_Z(-ERA), then 
 * dMdt = Q * dRdt *W
 *
 * @param[in] tt Epoch of request, in TT scale.
 * @param[in] eop An EopRecord instance, which will be used to extract: pole
 *            motion, i.e. xp, yp, corrections to CIP, i.e. dX and dY, 
 *            Delta UT1 and LOD.
 * @param[out] dMdt The derivative of M, where M is used to transform in the 
 *            sense r[TRS] = M * r[CRS], considering only Earth rotation. See
 *            also detail::dRdt.
 * @return The rotation quaternion to transform a vector in the sense:
 *         [TRS] = q * [CRS]
 *
 * @code
 *  // vectors r, v are in GCRF and we want the velocity vector in ITRF:
 *  Eigen::Matrix3d dMdt;
 *  const auto q_c2i = dso::c2i06a(mjd, eops, dMdt);
 *  Eigen::Vector3d v1 = q_c2i * v + dMdt * r;
 *  // inverse transformation, i.e. from ITRF to GCRF:
 *  Eigen::Vector3d v11 = q_c2i.conjugate() * v1 + dMdt.transpose() * r1;
 * @endcode
 */
Eigen::Quaterniond c2i06a(const MjdEpoch &tt, const EopRecord &eop,
                          Eigen::Matrix3d &dMdt) noexcept;

/** @brief Get the unit quaternion to transform from ITRS to GCRS, i.e.
 *  [TRS] = q * [CRS]
 *
 * The function will extract (observed) EOP data from the passed in eops
 * instance, i.e. (a) xp, (b) yp, (c) dut, amd (d) dXcip and dYcip. The 
 * remainind data needed to perform the computation (i.e. (X,Y)_CIP, CIO 
 * locator s, and TIO locator s') will be evaluated at function call, based 
 * on the passed in epoch (tt).
 *
 * @param[in] tt Epoch of request, in TT scale
 * @param[in] eops An EopRecord instance, holding EOP data for the epoch of
 *             request (i.e. tt). The following EOP data will be requested
 *             from the instance: (a) xp, (b) yp, and (c) dut.
 * @return A unit quaternion q, acting in the sense: [TRS] = q * [CRS]
 */
Eigen::Quaterniond c2i06a_bz(const MjdEpoch &tt, const EopRecord &eops) noexcept;

Eigen::Vector3d c2i06a(const MjdEpoch &tt, const EopRecord &eop, Eigen::Quaterniond &qc2tirs, Eigen::Quaterniond &qtirs2i) noexcept;

namespace detail {

/** @brief Trnaform (X,Y)_{CIP} coordinates to spherical angles E and d.
 *
 * X_cip, Y_cip are the the x, y coordinates of the CIP unit vector in the
 * GCRS. E and d are the coordinates of the CIP in the GCRS such that:
 * X = sin d cos E, Y = sin d sin E, Z = cos d
 * See IERS 2010, Eq. 5.7
 *
 * @param[in] xcip X coordinate of CIP, in [rad]
 * @param[in] ycip Y coordinate of CIP, in [rad]
 * @param[out] e Spherical angle E in [rad]
 * @param[out] d Spherical angle d in [rad]
 */
inline void xycip2spherical(double Xcip, double Ycip, double &d,
                            double &e) noexcept {
  /* Obtain the spherical angles E and d:
   * x = sin d cos E,
   * y = sin d sin E,
   * [Z = cos d]
   */
  const double r2 = Xcip * Xcip + Ycip * Ycip;
  e = (r2 > 0e0) ? std::atan2(Ycip, Xcip) : 0e0;
  d = std::atan(std::sqrt(r2 / (1e0 - r2)));
}

/** @brief Transformation matrix for the celestial motion of the CIP (relating
 * CIRS and GCRS).
 *
 * The expression for the transformation matrix for the celestial motion of the
 * CIP is taken from the IERS 2010 standards, as:
 * Q = R3 (−E) x R2 (−d) x R3 (E) x R3 (s)
 * Note that this function will return the rotation transformation C = Q^T
 *
 * The transformation arises from the motion of the CIP in the GCRS (i.e.
 * relating CIRS and GCRS), E and d being such that the coordinates of the CIP
 * in the GCRS are: x = sind cosE, y = sind sinE, z = cosd. x and y here, denote
 * the coordinates of the CIP in the GCRS to be used for the parameters (see
 * IERS 2010, Sec. 5.5.4) They can be computed from the function dso::xycip06a.
 *
 * s is a quantity, named “CIO locator”, which provides the position of the CIO
 * on the equator of the CIP corresponding to the kinematical deﬁnition of the
 * NRO in the GCRS when the CIP is moving with respect to the GCRS, between the
 * reference epoch and the date t due to precession and nutation. It can be
 * computed via a call to dso::s06
 *
 * The resulting quaternion gives the GCRS-to-CIRS transformation, i.e.
 * r[CIRS] = C * r[GCRS]
 *
 * @param[in] x X-CIP coordinate [rad]
 * @param[in] y Y-CIP coordinate [rad]
 * @param[in] s CIO in [rad]
 * @return Rotation quaternion q, to transform between CIRS and GCRS, in the
 * sense: r[CIRS] = q * r[GCRS]
 */
inline auto C(double Xcip, double Ycip, double s) noexcept {
  double d, e;
  xycip2spherical(Xcip, Ycip, d, e);
  /* C = [R3 (−E) x R2 (−d) x R3 (E) x R3 (s)]^T */
  return (Eigen::AngleAxisd(e + s, Eigen::Vector3d::UnitZ()) *
          Eigen::AngleAxisd(-d, Eigen::Vector3d::UnitY()) *
          Eigen::AngleAxisd(-e, Eigen::Vector3d::UnitZ()));
}

/** @brief Transformation matrix for the celestial motion of the CIP (relating
 * CIRS and GCRS).
 *
 * This version computes the rotation quaternion as:
 * C = [PN(X,Y) * R3(s)]^T
 * To compute this, we follow Eq.(7) and Eq. (A2) from Bizouard and Cheng 2023.
 *
 * The resulting quaternion gives the GCRS-to-CIRS transformation, i.e.
 * r[CIRS] = C * r[GCRS]
 *
 * See Bizouard, C., Cheng, Y. The use of the quaternions for describing the
 * Earth’s rotation. J Geod 97, 53 (2023).
 * https://doi.org/10.1007/s00190-023-01735-z
 *
 * @param[in] x X-CIP coordinate [rad]
 * @param[in] y Y-CIP coordinate [rad]
 * @param[in] s CIO in [rad]
 * @return Rotation quaternion q, to transform between CIRS and GCRS, in the
 * sense: r[CIRS] = q * r[GCRS]
 */
inline auto C_qimpl(double Xcip, double Ycip, double s) noexcept {
  double d, e;
  xycip2spherical(Xcip, Ycip, d, e);
  const double z = std::cos(d);
  const double cs2 = std::cos(s / 2e0);
  const double ss2 = std::sin(s / 2e0);
  /* components of the quaternion R3(s) * q_PN, where
   * 1. q_pn is extracted from Eq.(7)
   * 2. and multiplication follows Eq. (A2)
   */
  const double f = 1e0 / std::sqrt(2e0 * (1e0 + z));
  const double w = f * (cs2 * (1e0 + z));
  const double a = f * (cs2 * Ycip - ss2 * (-Xcip));
  const double b = f * (cs2 * (-Xcip) + ss2 * Ycip);
  const double c = f * (ss2 * (1e0 + z));
  return Eigen::Quaterniond(w, a, b, c);
}

/** @brief Transformation matrix for the celestial motion of the CIP ( relating
 * CIRS and GCRS).
 *
 * This version computes the rotation matrix as:
 * C = [PN(X,Y) * R3(s)]^T
 * To compute this, we follow Eq.(5.10) from IERS 2010.
 *
 * The resulting quaternion gives the GCRS-to-CIRS transformation, i.e.
 * r[CIRS] = C * r[GCRS]
 *
 * @note this version seems to provide a bit of less accuracy, see
 * test/sofa/rot_pn.cpp, case A.
 *
 * @param[in] x X-CIP coordinate [rad]
 * @param[in] y Y-CIP coordinate [rad]
 * @param[in] s CIO in [rad]
 * @return Rotation quaternion q, to transform between CIRS and GCRS, in the
 * sense: r[CIRS] = q * r[GCRS]
 */
inline auto C_rxyimpl(double x, double y, double s) noexcept {
  const double a = 1. / 2. + 1. / 8. * (x * x + y * y);
  Eigen::Matrix<double, 3, 3> PN;
  PN(0, 0) = 1e0 - a * x * x;
  PN(0, 1) = -a * x * y;
  PN(0, 2) = -x;
  PN(1, 0) = -a * x * y;
  PN(1, 1) = 1e0 - a * y * y;
  PN(1, 2) = -y;
  PN(2, 0) = x;
  PN(2, 1) = y;
  PN(2, 2) = 1e0 - a * (x * x + y * y);
  return Eigen::AngleAxisd(s, Eigen::Vector3d::UnitZ()) * PN;
}

/** @brief Compute the derivative of Earth rotation matrix as 
 * dRdt = d(R_z(-ERA))/dt.
 *
 * The (original) rotation R, transforms in the sense:
 * [CIRS] = R * [TIRS]
 *
 * Here we consider that d(ERA)/dt = omega. The matrix is thus computed as:
 *         | -omega * sin(ERA)   omega * cos(ERA)  0 |
 * dR/dt = | -omega * cos(ERA)  -omega * sin(ERA)  0 |
 *         |         0                 0           0 |
 *
 * @param[in] era Earth rotation angle (ERA) in [rad]
 * @return Earth rotation quaternion q, with q = R_z(-ERA).
 */
inline auto dRdt(double era, double omega) noexcept {
  const double os = omega * std::sin(era);
  const double oc = omega * std::cos(era);
  Eigen::Matrix3d dRdt = Eigen::Matrix3d::Zero();
  dRdt(0,0) = -os;
  dRdt(1,0) = -oc;
  dRdt(0,1) =  oc;
  dRdt(1,1) = -os;
  return dRdt;
}

/** @brief Compute the quaternion to describe the GCRS to ITRS transformation: 
 * [TRS] = q * [CRS]
 *
 * This implementation is based on the Bizouard and Cheng (2023) model. It
 * uses analytical expressions to compute all entries of a quaternion that
 * describes the rotation (transformation) between the ITRS and the GCRS,
 * given a series of input parameters.
 * See the formulae (13) and (14) in the provided reference.
 *
 * @param[in] era  Earth Rotation Angle [rad]
 * @param[in] s    The CIO locator [rad]
 * @param[in] sp   The TIO locator [rad] (i.e. s')
 * @param[in] Xcip X-component of the CIP in the GCRS, [rad]. For precise
 *                 applications, this value should be the one computed by the
 *                 IAU 2006/2000A model and corrected via the IERS published
 *                 values (i.e. EOP data). That is:
 *                 Xcip = X(IAU 2006/2000A) + δX
 * @param[in] Ycip Y-component of the CIP in the GCRS, [rad]. For precise
 *                 applications, this value should be the one computed by the
 *                 IAU 2006/2000A model and corrected via the IERS published
 *                 values (i.e. EOP data). That is:
 *                 Ycip = Y(IAU 2006/2000A) + δY
 * @param[in] xp   X-coordinate of the "polar coordinates", i.e. of the
 *                 Celestial Intermediate Pole (CIP) in the ITRS, [rad]
 * @param[in] yp   Y-coordinate of the "polar coordinates", i.e. of the
 *                 Celestial Intermediate Pole (CIP) in the ITRS, [rad]
 *
 * Bizouard, C., Cheng, Y. The use of the quaternions for describing the
 * Earth’s rotation. J Geod 97, 53 (2023).
 * https://doi.org/10.1007/s00190-023-01735-z
 */
Eigen::Quaterniond gcrs2itrs_quaternion(double era, double s, double sp,
                                        double Xcip, double Ycip, double xp,
                                        double yp) noexcept;

/** @brief Compute the GCRS to ITRS transformation (rotation) quaternion:
 * [TRS] = q * [CRS]
 *
 * This implementation is based on the IERS 2010 document, following the
 * so called 'CIO-based' transformation.
 *
 * @param[in] era  Earth Rotation Angle [rad]
 * @param[in] s    The CIO locator [rad]
 * @param[in] sp   The TIO locator [rad] (i.e. s')
 * @param[in] d    d-component (spherical angle) of the CIP in the GCRS, [rad]. 
 * @param[in] e    E-component (spherical angle) of the CIP in the GCRS, [rad]. 
 *                 For precise applications, this value should be the one 
 *                 computed by the IAU 2006/2000A model and corrected via the 
 *                 IERS published values (i.e. EOP data). That is:
 *                 Xcip = X(IAU 2006/2000A) + δX and
 *                 Ycip = Y(IAU 2006/2000A) + δY
 *                 To obtain spherical angles (E, d) from (X, Y)_{CIP}:
 *                 X = sin d cos E, Y = sin d sin E, Z = cos d
 *                 See function detail::xycip2spherical
 * @param[in] xp   X-coordinate of the "polar coordinates", i.e. of the
 *                 Celestial Intermediate Pole (CIP) in the ITRS, [rad]
 * @param[in] yp   Y-coordinate of the "polar coordinates", i.e. of the
 *                 Celestial Intermediate Pole (CIP) in the ITRS, [rad]
 * @return Rotation quaternion q, such that: [TRS] = q * [CRS]
 */
inline Eigen::Quaterniond c2i(double era, double s, double sp, double d,
                              double e, double xp, double yp) noexcept {
  using namespace Eigen;
  return AngleAxisd(yp, Vector3d::UnitX()) * AngleAxisd(xp, Vector3d::UnitY()) *
         AngleAxisd(-sp - era + e + s, Vector3d::UnitZ()) *
         AngleAxisd(-d, Vector3d::UnitY()) * AngleAxisd(-e, Vector3d::UnitZ());
}

/** @brief Compute the GCRS to TIRS transformation (rotation) quaternion:
 * [TIRS] = q * [CRS]
 *
 * This implementation is based on the IERS 2010 document, following the
 * so called 'CIO-based' transformation.
 *
 * @param[in] era  Earth Rotation Angle [rad]
 * @param[in] s    The CIO locator [rad]
 * @param[in] sp   The TIO locator [rad] (i.e. s')
 * @param[in] d    d-component (spherical angle) of the CIP in the GCRS, [rad]. 
 * @param[in] e    E-component (spherical angle) of the CIP in the GCRS, [rad]. 
 *                 For precise applications, this value should be the one 
 *                 computed by the IAU 2006/2000A model and corrected via the 
 *                 IERS published values (i.e. EOP data). That is:
 *                 Xcip = X(IAU 2006/2000A) + δX and
 *                 Ycip = Y(IAU 2006/2000A) + δY
 *                 To obtain spherical angles (E, d) from (X, Y)_{CIP}:
 *                 X = sin d cos E, Y = sin d sin E, Z = cos d
 *                 See function detail::xycip2spherical
 * @param[in] xp   X-coordinate of the "polar coordinates", i.e. of the
 *                 Celestial Intermediate Pole (CIP) in the ITRS, [rad]
 * @param[in] yp   Y-coordinate of the "polar coordinates", i.e. of the
 *                 Celestial Intermediate Pole (CIP) in the ITRS, [rad]
 * @return Rotation quaternion q, such that: [TRS] = q * [CRS]
 */
inline Eigen::Quaterniond c2tirs(double era, double s, double d,
                                 double e) noexcept {
  using namespace Eigen;
  return AngleAxisd(-era + e + s, Vector3d::UnitZ()) *
         AngleAxisd(-d, Vector3d::UnitY()) * AngleAxisd(-e, Vector3d::UnitZ());
}

/** @brief Compute the TIRS to ITRS transformation (rotation) quaternion:
 * [TRS] = q * [TIRS]
 *
 * This implementation is based on the IERS 2010 document, following the
 * so called 'CIO-based' transformation.
 *
 * @param[in] sp   The TIO locator [rad] (i.e. s')
 * @param[in] xp   X-coordinate of the "polar coordinates", i.e. of the
 *                 Celestial Intermediate Pole (CIP) in the ITRS, [rad]
 * @param[in] yp   Y-coordinate of the "polar coordinates", i.e. of the
 *                 Celestial Intermediate Pole (CIP) in the ITRS, [rad]
 * @return Rotation quaternion q, such that: [TRS] = q * [TIRS]
 */
inline Eigen::Quaterniond tirs2i(double xp, double yp, double sp) noexcept {
  using namespace Eigen;
  return AngleAxisd(yp, Vector3d::UnitX()) * AngleAxisd(xp, Vector3d::UnitY()) *
         AngleAxisd(-sp, Vector3d::UnitZ());
}

} /* namespace detail */

/** @brief Earth's roatation rate (omega, ω) in [rad/sec]
 *
 * @param[in] dlod Delta LOD (e.g. interpolated from an EOP series) in [sec]
 * @return Instantaneous Earth's roatation rate in [rad/sec]
 */
inline double earth_rotation_rate(double dlod) noexcept {
  return ::iers2010::OmegaEarth * (1e0 - dlod / 86400e0);
}

inline Eigen::Vector3d earth_rotation_axis(double lod) noexcept {
  return Eigen::Vector3d(0e0, 0e0, earth_rotation_rate(lod));
}

} /* namespace dso */

#endif
