/** @file
 * Declare a StokesCoefficients class, that can hold C and S coefficients for
 * a given degree and order, to assist handling of Spherical Harmonics.
 */

#ifndef __STOKES_HARMONIC_POTENTIAL_COEFFICIENTS_HPP__
#define __STOKES_HARMONIC_POTENTIAL_COEFFICIENTS_HPP__

#include "coeff_matrix_2d.hpp"
#include "iersconst.hpp"
#include <cassert>
#include "eigen3/Eigen/Eigen"

namespace dso {

/** A class to hold Stokes Coefficients C and S.
 *
 * Note that the Stokes coefficients ares stored in Lower Triangular but 
 * RECTANGULAR matrices, i.e. matrices that have the same number of rows and 
 * columns. Hence, desite the fact that each Stokes Coefficients instance 
 * has two members to designate the number of coefficients (m_degree and 
 * m_column), only one (m_rows) is relevant for storage. _Cnm and _Snm, 
 * know nothing about degree!
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

  /** Get the underlying data of the C coefficient matrix for given index.
   * Use with care! There is not certain way of telling the (n,m) index 
   * of the element requested. This function should only be used if we want 
   * to do smthng with all the elements of C, which is very rare.
   * This function is mstly used for assistance within the Expression 
   * Templates.
   */
  double Cdata(int i) const noexcept {return _Cnm.data()[i];}

  /** Get the underlying data of the S coefficient matrix for given index.
   * Use with care! There is not certain way of telling the (n,m) index 
   * of the element requested. This function should only be used if we want 
   * to do smthng with all the elements of S, which is very rare.
   * This function is mstly used for assistance within the Expression 
   * Templates.
   */
  double Sdata(int i) const noexcept {return _Snm.data()[i];}

public:
  /** Copy constructor */
  StokesCoeffs(const StokesCoeffs &other) noexcept
      : _GM(other._GM), _Re(other._Re), _cnormalized(other._cnormalized),
        m_degree(other.m_degree), m_order(other.m_order), _Cnm(other._Cnm),
        _Snm(other._Snm) 
  {
  }

  /** Move constructor */
  StokesCoeffs(StokesCoeffs &&other) noexcept
      : _GM(other._GM), _Re(other._Re), _cnormalized(other._cnormalized),
        m_degree(other.m_degree), m_order(other.m_order),
        _Cnm(std::move(other._Cnm)), _Snm(std::move(other._Snm)) 
  {
  }

  /** Assignment operator */
  StokesCoeffs &operator=(const StokesCoeffs &other) noexcept {
    if (this != &other) {
      _GM = other._GM;
      _Re = other._Re;
      _cnormalized = other._cnormalized;
      m_degree = other.m_degree;
      m_order = other.m_order;
      _Cnm = other._Cnm;
      _Snm = other._Snm;
    }
    return *this;
  }
  
  /** Move assignment operator */
  StokesCoeffs &operator=(StokesCoeffs &&other) noexcept {
    if (this != &other) {
      _GM = other._GM;
      _Re = other._Re;
      _cnormalized = other._cnormalized;
      m_degree = other.m_degree;
      m_order = other.m_order;
      _Cnm = std::move(other._Cnm);
      _Snm = std::move(other._Snm);
    }
    return *this;
  }

  /** Destructor is a no-op */
  ~StokesCoeffs() noexcept {};

  void swap(StokesCoeffs &b) noexcept {
    using std::swap;
    std::swap(_GM, b._GM);
    std::swap(_Re, b._Re);
    std::swap(_cnormalized, b._cnormalized);
    std::swap(m_degree, b.m_degree);
    std::swap(m_order, b.m_order);
    swap(_Cnm, b._Cnm);
    swap(_Snm, b._Snm);
  }

  /** Default constructor */
  StokesCoeffs() noexcept
      : _GM(::iers2010::GMe), _Re(::iers2010::Re), _cnormalized(true), m_degree(0),
        m_order(0), _Cnm{0, 0}, _Snm{0, 0} {}

  /** Constructor given degree, order, GM and radius R */
  StokesCoeffs(int n, int m=-1, double GM=::iers2010::GMe, double Re=::iers2010::Re)
      : _GM(GM), _Re(Re), _cnormalized(true), m_degree(n), 
        m_order((m > 0) ? m : n), _Cnm(n + 1, n + 1), _Snm(n + 1, n + 1) {}

  /* @brief Resize to new degree/order, deleting currrent data.
   *
   * Check current capacity and only re-allocated data if needed. 
   * 
   * @warning Note that data stored in the instance will be lost! If you need 
   *          to resize the instance, but keep the data already stored in it, 
   *          use the cresize function.
   */
  void resize(int degree, int order) noexcept;
  
  /* @brief Resize to new degree/order, keeping current data.
   *
   * Check current capacity and only re-allocated data if needed. 
   */
  void cresize(int degree, int order) noexcept;

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
      _Cnm.fill_with(0e0);
      _Snm.fill_with(0e0);
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
   * The degree and order of the passed in instance (\p sc), should be equal 
   * to or less than the ones of the calling instance.
   *
   * @warning calling instance should be larger than argument (\p sc).
   */
  StokesCoeffs &operator+=(const StokesCoeffs &sc);

  /** @brief Copy C and S coefficients in 'packed storage' to an Eigen::Matrix
   *
   * C will be recorded in the lower triangular part of the resulting (anm) 
   * matrix, including the main diagonal. S will be stored transposed in the 
   * upper triangular part of anm. Note that we assume that all Snm coeffs 
   * are equal to zero when m=0; hence these values are not stored.
   *
   * The returned matrix is of size (m_degree+1)x(m_degree+1).
   * 
   * Indexing between anm, C and S is as follows:
   * C(row, col) = anm(row, col), assuming row>=col
   * S(row, col) = anm(col, row) if row>cols and 0 if row==col, again 
   *               assuming row>=col
   *
   * @snippet ../test/unit_tests/stokes2anm.cpp example_usage
   */
  Eigen::MatrixXd to_anm() const noexcept {
    Eigen::MatrixXd anm = Eigen::MatrixXd::Zero(m_degree+1, m_degree+1);
    for (int c=0; c<=m_degree; c++) {
      for (int r=0; r<c; r++) {
        anm(r,c) = _Snm(c,r+1);
      }
      for (int r=c; r<=m_degree; r++) {
        anm(r,c) = _Cnm(r,c);
      }
    }
    return anm;
  }
  
  static StokesCoeffs from_anm(const Eigen::MatrixXd& anm, int max_order=-1) noexcept {
    if (max_order<0) max_order = anm.rows()-1;
    int max_degree = anm.rows()-1;
    StokesCoeffs cs(max_degree, max_order);
    for (int m = 0; m <= max_degree; m++) {
        for (int l=0; l<m; l++) {
          cs.S(m,l+1) = anm(l,m);
        // (1,1) -> (0,1) m = 1
        // (2,1) -> (0,2) m = 2
        // (2,2) -> (1,2) m = 2
        // (3,1) -> (0,3) m = 3
        // (3,2) -> (1,3) m = 3
        // (3,3) -> (2,3) m = 3
        }
        for (int l = m; l <= max_degree; l++) {
          cs.C(l, m) = anm(l, m);
        }
    }
    return cs;
  }

}; /* StokesCoeffs */

inline void swap(StokesCoeffs &a, StokesCoeffs &b) noexcept {
  a.swap(b);
}

} /* namespace dso */

#endif
