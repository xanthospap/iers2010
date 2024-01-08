/** @file
 * Declare a StokesCoefficients class, that can hold C and S coefficients for
 * a given degree and order, to assist handling of Spherical Harmonics.
 */

#ifndef __STOKES_HARMONIC_POTENTIAL_COEFFICIENTS_HPP__
#define __STOKES_HARMONIC_POTENTIAL_COEFFICIENTS_HPP__

#include "coeff_matrix_2d.hpp"
#include "iersconst.hpp"
#include <cassert>

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

  /** Expression Template: Assist the summation of two StokesCoeffs instances 
   * using a proxy object. This proxy only holds references to the underlying 
   * StokesCoeffs instances, using lazy evaluation at the time od assignment.
   *
   * _SumProxy instances can only handle StokesCoeffs instances of same size 
   * (degree and order) and parameters (i.e. GM, Re, and normalized).
   */
  template<typename T1, typename T2>
  struct _SumProxy {
    const T1 &lhs;
    const T2 &rhs;
    int max_degree() const noexcept {return lhs.max_degree();}
    int max_order() const noexcept {return lhs.max_order();}
    double GM() const noexcept {return lhs.GM();}
    double Re() const noexcept {return lhs.Re();}
    bool normalized() const noexcept {return lhs.normalized();}

    double C(int i, int j) const noexcept {
      return lhs.C(i,j) + rhs.C(i,j);
    }
    double S(int i, int j) const noexcept {
      return lhs.S(i,j) + rhs.S(i,j);
    }
    double Cdata(int i) const noexcept {return lhs.Cdata(i) + rhs.Cdata(i);}
    double Sdata(int i) const noexcept {return lhs.Sdata(i) + rhs.Sdata(i);}
    
    _SumProxy(const T1 &ml, const T2 &mr) noexcept : lhs(ml), rhs(mr) {
      assert((ml.max_degree() == mr.max_degree()) &&
             (ml.max_order() == mr.max_order()) && (ml.GM() == mr.GM()) &&
             (ml.Re() == mr.Re()) && (ml.normalized() == mr.normalized()));
    }

    /** Add two _SumProxy instances */
    template <typename U1, typename U2>
    _SumProxy<_SumProxy, _SumProxy<U1, U2>>
    operator+(const _SumProxy<U1, U2> &other) const noexcept {
      return _SumProxy<_SumProxy, _SumProxy<U1, U2>>(*this, other);
    }

  }; /* _SumProxy */

  /** Expression Template: Assist the scaling (i.e. multiplication) of a 
   * StokesCoeffs instance using a proxy object. This proxy only holds a
   * reference to the underlying StokesCoeffs instance and a scaling factor 
   * (double), using lazy evaluation at the time od assignment.
   */
  template<typename T>
  struct _ScaledProxy {
    const T &lhs;
    double scale;
    int max_degree() const noexcept {return lhs.m_degree;}
    int max_order() const noexcept {return lhs.m_order;}
    double GM() const noexcept {return lhs._GM;}
    double Re() const noexcept {return lhs._Re;}
    bool normalized() const noexcept {return lhs._cnormalized;}
    double C(int i, int j) const noexcept {
      return lhs.C(i,j) * scale;
    }
    double S(int i, int j) const noexcept {
      return lhs.S(i,j) * scale;
    }
    double Cdata(int i) const noexcept {return lhs.Cdata(i)*scale;}
    double Sdata(int i) const noexcept {return lhs.Sdata(i)*scale;}

    _ScaledProxy(const T&m, double f) noexcept : lhs(m), scale(f) {}

    /** Add two _ScaledProxy instances to create a _SumProxy */
    _SumProxy<_ScaledProxy, _ScaledProxy>
    operator+(const _ScaledProxy &other) const noexcept {
      return _SumProxy<_ScaledProxy, _ScaledProxy>(*this, other);
    }

  }; /* _ScaledProxy */

public:
  /** Copy constructor */
  StokesCoeffs(const StokesCoeffs &other) noexcept
      : _GM(other._GM), _Re(other._Re), _cnormalized(other._cnormalized),
        m_degree(other.m_degree), m_order(other.m_order), _Cnm(other._Cnm),
        _Snm(other._Snm) 
  {
  }

  template<typename T> StokesCoeffs(T &&other) noexcept 
      : _GM(other.GM()), _Re(other.Re()), _cnormalized(other.normalized()),
        m_degree(other.max_degree()), m_order(other.max_order()), _Cnm(m_degree+1, m_order+1),
        _Snm(m_degree+1, m_order+1) 
  {
    const int N = _Cnm.num_elements();
    double *__restrict__ c = _Cnm.data();
    for (int i=0; i<N; i++) 
      c[i] = other.Cdata(i);
    double *__restrict__ s = _Snm.data();
    for (int i=0; i<N; i++) 
      s[i] = other.Sdata(i);
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
  //StokesCoeffs() noexcept
  //    : _GM(iers2010::GMe), _Re(iers2010::Re), _cnormalized(true), m_degree(0),
  //      m_order(0), _Cnm{0, 0}, _Snm{0, 0} {}

  /** Constructor given degree, order, GM and radius R */
  StokesCoeffs(int n, int m, double GM, double Re)
      : _GM(GM), _Re(Re), _cnormalized(true), m_degree(n), m_order(m),
        _Cnm(n + 1, m + 1), _Snm(n + 1, m + 1) 
  {
  }

  /** Constructor given degree (n)*/
  //StokesCoeffs(int n)
  //    : _GM(iers2010::GMe), _Re(iers2010::Re), _cnormalized(true), m_degree(n),
  //      m_order(n), _Cnm(n + 1, n + 1), _Snm(n + 1, n + 1) {}

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
   * The degree and order of the passed in instance (\p sc), should be equal 
   * to or less than the ones of the calling instance.
   *
   * @warning calling instance should be larger than argument (\p sc).
   */
  StokesCoeffs &operator+=(const StokesCoeffs &sc);

  template<typename T>
  _SumProxy<StokesCoeffs, T> operator+(const T &rhs) const noexcept {
    return _SumProxy<StokesCoeffs, T>(*this, rhs);
  }

  /** Overload += operator. The right hand side can be a _SumProxy or a
   * _ScaledProxy.
   *
   * Note that the right hand side can either:
   * 1. have the same size (i.e. max degree and order) with the calling 
   *    instance, or
   * 2. have smaller size (i.e. lhs.degree >= rhs.degree and 
   *    lhs.order >= lhs.order0
   */
  template <typename T> StokesCoeffs &operator+=(const T &rhs) noexcept {
    assert((this->max_degree() >= rhs.max_degree()) &&
           (this->max_order() >= rhs.max_order()) && (this->GM() == rhs.GM()) &&
           (this->Re() == rhs.Re()) &&
           (this->normalized() == rhs.normalized()));
    if ((this->max_degree() == rhs.max_degree()) &&
        (this->max_order() == rhs.max_order())) {
      const int N = _Cnm.num_elements();
      double *__restrict__ c = _Cnm.data();
      for (int i = 0; i < N; i++)
        c[i] += rhs.Cdata(i);
      double *__restrict__ s = _Snm.data();
      for (int i = 0; i < N; i++)
        s[i] += rhs.Sdata(i);
    } else {
      for (int m = 0; m <= rhs.max_order(); m++) {
        for (int n = m; n <= rhs.max_degree(); n++) {
          C(n, m) += rhs.C(n, m);
          S(n, m) += rhs.S(n, m);
        }
      }
    }
    return *this;
  }

}; /* StokesCoeffs */

inline void swap(StokesCoeffs &a, StokesCoeffs &b) noexcept {
  a.swap(b);
}

/** _SumProxy <- StokesCoeffs + _SumProxy */
template <typename T1, typename T2>
inline StokesCoeffs::_SumProxy<T1, T2> operator+(const T1 &lhs,
                                                 const T2 &rhs) noexcept {
  return StokesCoeffs::_SumProxy<T1, T2>(lhs, rhs);
}

/** _ScaledProxy <- s * StokesCoeffs */
inline StokesCoeffs::_ScaledProxy<StokesCoeffs>
operator*(double f, const StokesCoeffs &sc) noexcept {
  return StokesCoeffs::_ScaledProxy<StokesCoeffs>(sc, f);
}

} /* namespace dso */

#endif
