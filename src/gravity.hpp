/** @file
 * Computation of gravity acceleration using a geopotential model.
 */

#ifndef __DSO_GEOPOTENTIAL_ACCELERATAION_HPP__
#define __DSO_GEOPOTENTIAL_ACCELERATAION_HPP__

#include "stokes_coefficients.hpp"
#include "eigen3/Eigen/Eigen"

namespace dso {

/** Acceleration due to point mass at r_cb on a mass at r.
 *
 * Compute the aceleration induced on a body at the position vector r, caused
 * by a point mass with gravitational constant GM at position rcb.
 * This is normally used to compute third-body acceleration on a satellite, 
 * induced by e.g. the Sun or Moon.
 *
 * @param[in] r Geocentric position of attracted mass, i.e. satellite in [m].
 *              The vector holds Cartesian components, i.e. r=(x,y,z).
 * @param[in] rcb Geocentric position vector of the attracting body, i.e. the 
 *              'third body' (e.g. Sun or Moon) in [m].
 * @param[in] GMcb gravitational constant of the 'third body', i.e. G * M_cb, 
 *              in [m^3/ sec^2]
 * @return Acceleration induced on the mass, in Cartesian components, in 
 *         units of [m/s^2]. I.e. a = (a_x, a_y, a_z)
 *
 * @note Instead of using [m] as input units, [km] can also be used, as long
 *       as it is used consistently for ALL inputs (r, rcb and GM). in this 
 *       case, the resulting acceleration will be given in units of [km/s^2].
 */
Eigen::Matrix<double, 3, 1>
point_mass_acceleration(const Eigen::Matrix<double, 3, 1> &r,
                        const Eigen::Matrix<double, 3, 1> &rcb,
                        double GMcb) noexcept;

/** Acceleration due to point mass at r_cb on a mass at r.
 *
 * Same as above, only in this case we also compute and return the Jacobian 
 * matrix da/dr, i.e.
 *
 *     | dax/dx dax/dy dax/dz |
 * J = | day/dx day/dy day/dz |
 *     | daz/dx daz/dy daz/dz |
 *
 * @param[in] r Geocentric position of attracted mass, i.e. satellite in [m].
 *              The vector holds Cartesian components, i.e. r=(x,y,z).
 * @param[in] rcb Geocentric position vector of the attracting body, i.e. the 
 *              'third body' (e.g. Sun or Moon) in [m].
 * @param[in] GMcb gravitational constant of the 'third body', i.e. G * M_cb, 
 *              in [m^3/ sec^2]
 * @param[out] jacobian The Jacobian 3x3 matrix da/dr
 * @return Acceleration induced on the mass, in Cartesian components, in 
 *         units of [m/s^2]. I.e. a = (a_x, a_y, a_z)
 *
 * @note Instead of using [m] as input units, [km] can also be used, as long
 *       as it is used consistently for ALL inputs (r, rcb and GM). in this 
 *       case, the resulting acceleration will be given in units of [km/s^2].
 */
Eigen::Matrix<double, 3, 1>
point_mass_acceleration(const Eigen::Matrix<double, 3, 1> &r,
                        const Eigen::Matrix<double, 3, 1> &rcb,
                        double GMcb, 
                        Eigen::Matrix<double,3,3> &jacobian) noexcept;

/** @brief Compute normalised associated Legendre functions Pnm
 *
 * The algorithm employed here to perform the computations is the 
 * "forward recursions" method, see Holmes et al, 2002.
 *
 * @param[in] theta Geocentric latitude for expansion in [rad].
 * @param[in] max_degree Max degree for expansion. i.e. n in Pnm, inclusive
 * @param[in] max_order  Max order for expansion, i.e. m in Pnm, inclusive
 * @param[out] Pnm       Matrix to store the computed Pnm values; must be 
 *                       at least large enough to hold Pnm's for 
 *                       n=[0,max_degree] and m=[0,max_order].
 * @return Always zero.
 *
 * Holmes, S., Featherstone, W. A unified approach to the Clenshaw summation 
 * and the recursive computation of very high degree and order normalised 
 * associated Legendre functions. Journal of Geodesy 76, 279–299 (2002). 
 * https://doi.org/10.1007/s00190-002-0216-2
 */
int normalised_associated_legendre_functions(
    double theta, int max_degree, int max_order,
    CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &Pnm) noexcept;

/** Spherical harmonics of Earth's gravity potential to acceleration and
 *  gradient using the algorithm due to Cunningham. The acceleration and
 *  gradient are computed in Cartesian components, i.e.
 *
 *  acceleration = (ax, ay, az), and
 *
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
