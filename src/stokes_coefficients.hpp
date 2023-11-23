/** @file
 * Declare a StokesCoefficients class, that can hold C and S coefficients for
 * a given degree and order, to assist handling of Spherical Harmonics.
 */

#ifndef __STOKES_HARMONIC_POTENTIAL_COEFFICIENTS_HPP__
#define __STOKES_HARMONIC_POTENTIAL_COEFFICIENTS_HPP__

#include "coeff_matrix_2d.hpp"
#include "iersconst.hpp"

namespace dso {

/** A class to hold Stokes Coefficients C and S.
 *
 * C and S coefficients are stored in Lower-Triangular matrices, stored in
 * Column-Major manner.
 * Except for Stokes Coefficients, the class also holds the gravitational
 * constant GM and the radius Re, to be used when (and if) a spherical
 * harmonics expansion is needed.
 *
 * For given degree (N) and order (M), the instance will hols the Stokes
 * coefficients Cnm for n in range [0, N].
 */
class StokesCoeffs {
private:
  /* gravitational constant times mass of Earth [m^3 s^-2] */
  double _GM;
  /* reference radius of the spherical harmonics [m] */
  double _Re;
  /* coefficients are normaliized (?) */
  bool _cnormalized;
  /* maximum degree */
  int m_degree;
  /* maximum order */
  int m_order;
  /* Cnm coefficients */
  CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> _Cnm;
  /* Snm coefficients */
  CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> _Snm;

public:
  /** Default constructor */
  StokesCoeffs() noexcept
      : _GM(iers2010::GMe), _Re(iers2010::Re), _cnormalized(true), m_degree(0),
        m_order(0), _Cnm{0, 0}, _Snm{0, 0} {}

  /** Constructor given degree, order, GM and radius R */
  StokesCoeffs(int n, int m, double GM, double Re)
      : _GM(GM), _Re(Re), _cnormalized(true), m_degree(n), m_order(m),
        _Cnm(n + 1, m + 1), _Snm(n + 1, m + 1) {}

  /** Constructor given degree (n)*/
  StokesCoeffs(int n)
      : _GM(iers2010::GMe), _Re(iers2010::Re), _cnormalized(true), m_degree(n),
        m_order(n), _Cnm(n + 1, n + 1), _Snm(n + 1, n + 1) {}

  /* @brief Resize; check current capacity and only re-allocated data if
   *      needed. m_degree set to new value.
   */
  void resize(int degree, int order) noexcept;

  /** shrink dimensions (i.e. degree and order) without touching capacity */
  int shrink_dimensions(int new_degree, int new_order) noexcept {
    if ((new_order <= new_degree) && (new_degree <= m_degree) &&
        (new_order <= m_order)) {
      m_degree = new_degree;
      m_order = new_order;
      return 0;
    }
    return 1;
  }

  /** get max degree */
  int max_degree() const noexcept { return m_degree; }

  /** get max order */
  int max_order() const noexcept { return m_order; }

  /** get gravitational constant times mass [kg^3/m^2] */
  double GM() const noexcept { return _GM; }

  /** get radious [m] */
  double Re() const noexcept { return _Re; }

  /** return true if Stokes coefficients are normalized */
  bool normalized() const noexcept { return _cnormalized; }

  /** get/set gravitational constant times mass [kg^3/m^2] */
  double &GM() noexcept { return _GM; }

  /** get/set radious [m] */
  double &Re() noexcept { return _Re; }

  /** get/set true/false depending on wether the Stokes coefficients are
   * normalized
   */
  bool &normalized() noexcept { return _cnormalized; }

  /** get the J2 term, i.e. -C_nm for n=2 and m=0 */
  double J2() const noexcept { return -_Cnm(2, 0); };

  /** set the Cnm and Snm coefficients to zero */
  void clear() noexcept {
    if (m_degree) {
      _Cnm.fill_with(0e0);
      _Snm.fill_with(0e0);
    }
  }

  /** scale the Cnm and Snm coefficients */
  void scale(double factor) noexcept {
    if (m_degree) {
      _Cnm.multiply(factor);
      _Snm.multiply(factor);
    }
  }

  /** get the Cnm coefficient */
  double C(int n, int m) const noexcept { return _Cnm(n, m); }

  /** get/set the Cnm coefficient */
  double &C(int n, int m) noexcept { return _Cnm(n, m); }

  /** get the Snm coefficient */
  double S(int n, int m) const noexcept { return _Snm(n, m); }

  /** get/set the Snm coefficient */
  double &S(int n, int m) noexcept { return _Snm(n, m); }

  /** get the Cnm coefficient matrix */
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &Cnm() noexcept {
    return _Cnm;
  }

  /** get the Cnm coefficient matrix */
  const CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &
  Cnm() const noexcept {
    return _Cnm;
  }

  /** get the Snm coefficient matrix */
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &Snm() noexcept {
    return _Snm;
  }

  /** get the Snm coefficient matrix */
  const CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &
  Snm() const noexcept {
    return _Snm;
  }

  /** @brief Add Stokes coefficients, not neccesserily of same size.
   *
   * The function will add the coefficients of the passed in instance (\p sc)
   * to the corresponding of the calling instance.
   * The degree and order of the passed in instance, should be equal to or
   * less than the ones of the calling instance.
   *
   * @warning calling instance should be larger than argument (\p sc)
   */
  StokesCoeffs &operator+=(const StokesCoeffs &sc);

}; /* StokesCoeffs */

} /* namespace dso */

#endif