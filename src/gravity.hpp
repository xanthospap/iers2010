/** @file
 * Computation of gravity acceleration using a geopotential model.
 */

#ifndef __DSO_GEOPOTENTIAL_ACCELERATAION_HPP__
#define __DSO_GEOPOTENTIAL_ACCELERATAION_HPP__

#include "stokes_coefficients.hpp"
#include "eigen3/Eigen/Eigen"

namespace dso {

/** Spherical harmonics of Earth's gravity potential to acceleration and
 *  gradient using the algorithm due to Cunningham. The acceleration and
 *  gradient are computed in Cartesian components, i.e.
 *  acceleration = (ax, ay, az), and
 *             | dax/dx dax/dy dax/dz |
 *  gradient = | day/dx day/dy day/dz |
 *             | daz/dx daz/dy daz/dz |
 *
 * @param[in] cs Normalized Stokes coefficients of spherical harmonics
 * @param[in] r  Position vector of satellite (aka point of computation) in
 *               an ECEF frame (e.g. ITRF)
 * @param[in] max_degree Max degree of spherical harmonics expansion. If set 
 *               to a negative number, the degree of expansion will be derived 
 *               from the cs input parameter, i.e. cs.max_degree()
 * @param[in] max_order  Max order of spherical harmonics expansion. If set 
 *               to a negative number, the order of expansion will be derived 
 *               from the cs input parameter, i.e. cs.max_order()
 * @param[in] Re Equatorial radius of the Earth in [m]. If set to a negative 
 *               number, then the cs.Re() method will be used to get it.
 * @param[in] GM Gravitational constant of Earth. If set to a negative 
 *               number, then the cs.GM() method will be used to get it.
 * @param[out] acc Acceleration in cartesian components in [m/s^2]
 * @param[out] gradient Gradient of acceleration in cartesian components
 * @param[in] W   A convinient storage space, as Column-Wise Lower Triangular
 *                matrix off dimensions at least (max_degree+2, max_degree+2).
 *                If not given, then the function will allocate and free the 
 *                required memmory.
 * @param[in] M   A convinient storage space, as Column-Wise Lower Triangular
 *                matrix off dimensions at least (max_degree+2, max_degree+2)
 *                If not given, then the function will allocate and free the 
 *                required memmory.
 * @return        Anything other than zero denotes an error.
 */
int sh2gradient_cunningham(
    const dso::StokesCoeffs &cs, const Eigen::Matrix<double, 3, 1> &r,
    Eigen::Matrix<double, 3, 1> &acc, Eigen::Matrix<double, 3, 3> &gradient,
    int max_degree=-1, int max_order=-1, double Re=-1, double GM=-1,
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> *W =
        nullptr,
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> *M =
        nullptr) noexcept;

} /* namespace dso */

#endif
