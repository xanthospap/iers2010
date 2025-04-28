/** @file
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

/** Get the unit quaternion to transform from ITRS to GCRS, i.e.
 *  r_{GCRS} = q * r_{ITRS}
 *
 * @param[in] tt Epoch of request, in TT scale
 * @param[in] eops An EopRecord instance, holding EOP data for the epoch of
 *               request (i.e. tt)
 * @param[out] fargs If not NULL, at output it will hold the luni-solar and
 *             planetary arguments used in the computation of (X,Y). Since we
 *             are computing them, we might as well return them! If not NULL,
 *             the array should be large enough to hold 14 doubles, i.e.
 *         [l, l', F, D, Om, L_Me, L_Ve, L_E, L_Ma, L_J, L_Sa, L_U, L_Ne, p_A]
 *             all in units of [rad].
 *
 * @return A unit quaternion q, acting in the sense: r_{GCRS} = q * r_{ITRS}
 */
Eigen::Quaterniond itrs2gcrs_quaternion(const MjdEpoch &tt,
                                        const EopRecord &eops,
                                        double *fargs = nullptr) noexcept;

/** Earth's roatation rate (omega, ω) in [rad/sec]
 *
 * @param[in] dlod Delta LOD (e.g. interpolated from an EOP series) in [sec]
 * @return Instantaneous Earth's roatation rate in [rad/sec]
 */
inline double earth_rotation_rate(double dlod) noexcept {
  return ::iers2010::OmegaEarth * (1e0 - dlod / 86400e0);
}

namespace detail {

/** @brief polar motion matrix W = [R3(-s') R2(xp) R1(yp)]^T
 *
 * The matrix operates in the sense: r_ITRS = W * r_TIRS
 * It is the transpose of the matrix described in Ch. 5.4.1 of the IERS 2010
 * standards.
 */
inline auto W(double xp, double yp, double sp) noexcept {
  /* W = [R3(-s') R2(xp) R1(yp)]^T */
  using namespace Eigen;
  return (AngleAxisd(yp, Vector3d::UnitX()) *
          AngleAxisd(xp, Vector3d::UnitY()) *
          AngleAxisd(-sp, Vector3d::UnitZ()));
}

inline auto C(double x, double y, double s) noexcept {
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

inline auto R(double era) noexcept {
  return Eigen::AngleAxisd(-era, Eigen::Vector3d::UnitZ());
}

/** @brief Compute the quaternion to describe the GCRS to ITRS transformation.
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

/** @brief Compute the GCRS to ITRS transformation (rotation) matrix.
 *
 * This implementation is based on the IERS 2010 document, following the
 * so called 'CIO-based' transformation.
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
 */
Eigen::Matrix<double, 3, 3> gcrs2itrs(double era, double s, double sp,
                                      double Xcip, double Ycip, double xp,
                                      double yp) noexcept;

/** @brief Compute the GCRS to ITRS transformation (rotation) matrix.
 */
Eigen::Matrix<double, 3, 3>
gcrs2itrs(double era, double s, double sp, double Xcip, double Ycip, double xp,
          double yp, double lod, Eigen::Matrix<double, 3, 3> &dRdt) noexcept;

/**
 */
Eigen::Matrix<double, 3, 3>
gcrs2tirs(double era, double s, double sp, double Xcip, double Ycip, double xp,
          double yp, Eigen::Matrix<double, 3, 3> *W = nullptr) noexcept;
} /* namespace detail */

/** @brief Return the polar motion rotation matrix W in the sense: TIRS-to-ITRS.
 *
 * Polar motion rotation matrix is computed by the formula:
 *  W = [R3(-s') R2(xp) R1(yp)]^(-1)
 * xp and yp being the “polar coordinates” of the Celestial Intermediate Pole
 * (CIP) in the ITRS and s′ being a quantity, named “TIO locator”, which
 * provides the position of the TIO on the equator of the CIP corresponding to
 * the kinematical deﬁnition of the “non-rotating” origin (NRO) in the ITRS when
 * the CIP is moving with respect to the ITRS due to polar motion.
 *
 * This matrix works as follows:
 * r_ITRS = W * r_TIRS
 *
 * @warning This matrix is the inverse/transpose of the matrix described in
 * Ch. 5.4.1 of the IERS 2010 standards.
 *
 * @param[in] xp Polar coordinate x in [rad]
 * @param[in] yp Polar coordinate y in [rad]
 * @param[in] sp TIO locator (s') in [rad]
 */
inline auto polar_motion_matrix(double xp, double yp, double sp) noexcept {
  return detail::W(xp, yp, sp);
}

} /* namespace dso */

#endif
