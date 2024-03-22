/** @file
 * Naive implementation of 2-Dimensional matrices; note that this
 * implementation targets 2-dimensional matrices that are only meant to store
 * data NOT perform arithmetic operations.
 */

#ifndef __COMPACT_2D_SIMPLE_MATRIX_HPP__
#define __COMPACT_2D_SIMPLE_MATRIX_HPP__

#include "coeff_matrix_storage.hpp"
#include <cassert>
#include <cstring>
#ifdef DEBUG
#include <cstdio>
#endif

namespace dso {

/** @brief A naive implementation of a 2-d dense matrix.
 *  The main objective of this class is to store data (e.g. coefficients) and
 *  not perform arithmetic operations.
 */
template <MatrixStorageType S> class CoeffMatrix2D {
private:
  StorageImplementation<S> m_storage; /** storage type; dictates indexing */
  double *m_data{nullptr};            /** the actual data */
  std::size_t _capacity{0}; /** number of doubles in allocated memory arena */

  /** Access an element from the underlying data; use with care IF needed */
  double data(int i) const noexcept { return m_data[i]; }

public:
  template <typename T1, typename T2> struct _SumProxy;

  /** Expression Template: Structure to hold a scaled CoeffMatrix2D (i.e. the
   * multiplication of some a matrix by a real number.
   */
  template <typename T1> struct _ScaledProxy {
    const T1 &mat;
    double fac;
    int rows() const noexcept { return mat.rows(); }
    int cols() const noexcept { return mat.cols(); }
    _ScaledProxy(const T1 &t1, double d) noexcept : mat(t1), fac(d) 
    {
      printf("[DEBUG] Called _ScaledProxy c'tor\n");
    };
    double operator()(int i, int j) const noexcept { return mat(i, j) * fac; }
    double data(int i) const noexcept { return mat.data(i) * fac; }
    template <typename T2>
    _SumProxy<_ScaledProxy, _ScaledProxy<T2>>
    operator+(const _ScaledProxy<T2> &s2) const noexcept {
      return _SumProxy<_ScaledProxy, _ScaledProxy<T2>>(*this, s2);
    }
  }; /* _ScaledProxy */

  template <typename T1, typename T2> struct _SumProxy {
    const T1 &lhs;
    const T2 &rhs;
    int rows() const noexcept { return lhs.rows(); }
    int cols() const noexcept { return lhs.cols(); }
    double operator()(int i, int j) const noexcept {
      return rhs(i, j) + lhs(i, j);
    }
    double data(int i) const noexcept { return lhs.data(i) + rhs.data(i); }
    _SumProxy(const T1 &t1, const T2 &t2) noexcept : lhs(t1), rhs(t2) {
      assert((lhs.rows() == rhs.rows()) && (lhs.cols() == rhs.cols()));
      printf("[DEBUG] Called _SumProxy c'tor\n");
    }
    /** Allow for _SumProxy<T1,T2> + _ScaledProxy<U> */
    template <typename U>
    _SumProxy<_SumProxy<T1, T2>, _ScaledProxy<U>>
    operator+(const _ScaledProxy<U> &scaled) const noexcept {
      return _SumProxy<_SumProxy<T1, T2>, _ScaledProxy<U>>(*this, scaled);
    }
    /** Allow for _SumProxy<T1,T2> + CoeffMatrix2D<S> */
    _SumProxy<_SumProxy<T1, T2>, CoeffMatrix2D>
    operator+(const CoeffMatrix2D &mat) const noexcept {
      return _SumProxy<_SumProxy<T1, T2>, CoeffMatrix2D>(*this, mat);
    }
  }; /* SumProxy */

  /** Swap current instance with another */
  void swap(CoeffMatrix2D<S> &b) noexcept {
    using std::swap;
    StorageImplementation<S> tmp_s = b.m_storage;
    std::size_t tmp_c = b._capacity;

    swap(m_data, b.m_data);
    b.m_storage = this->m_storage;
    b._capacity = this->_capacity;
    this->m_storage = tmp_s;
    this->_capacity = tmp_c;

    return;
  }

  /** get number of rows */
  constexpr int rows() const noexcept { return m_storage.nrows(); }

  /** get number of columns */
  constexpr int cols() const noexcept { return m_storage.ncols(); }

  /** get number of elements */
  constexpr std::size_t num_elements() const noexcept {
    return m_storage.num_elements();
  }

  /** @brief Element indexing (rows and columns start from 0 --not 1--) */
  double &operator()(int i, int j) noexcept {
#ifdef DEBUG
    assert(i < rows() && j < cols());
    assert(m_storage.index(i, j) >= 0 &&
           m_storage.index(i, j) < (int)m_storage.num_elements());
#endif
    return m_data[m_storage.index(i, j)];
  }

  /** @brief Element indexing (rows and columns start from 0 --not 1--) */
  const double &operator()(int i, int j) const noexcept {
#ifdef DEBUG
    assert(i < rows() && j < cols());
    assert(m_storage.index(i, j) >= 0 &&
           m_storage.index(i, j) < (int)m_storage.num_elements());
#endif
    return m_data[m_storage.index(i, j)];
  }

  /** @brief Row/Column indexing (rows and columns start from 0 --not 1--)
   *
   * If the data is stored in a Row-Wise manner, this function will return
   * the offset of the ith row; if data is stored in a column-wise manner,
   * it will return the index of the first elelement of the ith column.
   */
  const double *slice(int i) const noexcept {
    return m_data + m_storage.slice(i);
  }

  /** @brief Row/Column indexing (rows and columns start from 0 --not 1--)
   *
   * If the data is stored in a Row-Wise manner, this function will return
   * the offset of the ith row; if data is stored in a column-wise manner,
   * it will return the index of the first elelement of the ith column.
   */
  double *slice(int i) noexcept { return m_data + m_storage.slice(i); }

  /** @brief Pointer to the begining of a given column.
   *
   * Only defined if the MatrixStorageType uses some kind of Column Major
   * storage sequence.
   */
  template <MatrixStorageType t = S,
            std::enable_if_t<StorageImplementation<t>::isColMajor, bool> = true>
  double *column(int j) noexcept {
    return slice(j);
  }

  /** @brief Const pointer to the begining of a given column.
   *
   * Only defined if the MatrixStorageType uses some kind of Column Major
   * storage sequence.
   */
  template <MatrixStorageType t = S,
            std::enable_if_t<StorageImplementation<t>::isColMajor, bool> = true>
  const double *column(int j) const noexcept {
    return slice(j);
  }

  /** @brief Pointer to the begining of a given row.
   *
   * Only defined if the MatrixStorageType uses some kind of Row Major
   * storage sequence.
   */
  template <MatrixStorageType t = S,
            std::enable_if_t<StorageImplementation<t>::isRowMajor, bool> = true>
  double *row(int j) noexcept {
    return slice(j);
  }

  /** @brief Const pointer to the begining of a given row.
   *
   * Only defined if the MatrixStorageType uses some kind of Row Major
   * storage sequence.
   */
  template <MatrixStorageType t = S,
            std::enable_if_t<StorageImplementation<t>::isRowMajor, bool> = true>
  const double *row(int j) const noexcept {
    return slice(j);
  }

  /** set all elements of the matrix equal to some value */
  void fill_with(double val) noexcept {
    std::fill(m_data, m_data + m_storage.num_elements(), val);
  }

  /** multiply all elements of matrix with given value */
  void multiply(double value) noexcept {
    std::transform(m_data, m_data + m_storage.num_elements(), m_data,
                   [=](double d) { return d * value; });
  }

  /** get a pointer to the data */
  const double *data() const noexcept { return m_data; }

  /** get a non-const pointer to the data */
  double *data() noexcept { return m_data; }

  /** Constructor using number of rows and columns; for some
   * MatrixStorageType's, the number of columns may not be needed.
   */
  CoeffMatrix2D(int rows, int cols) noexcept
      : m_storage(rows, cols), m_data(new double[m_storage.num_elements()]),
        _capacity(m_storage.num_elements()) {
    printf("[DEBUG] CoeffMatrix2D c'tor; mem at: %p size=%ld\n", (void*)m_data, _capacity);
#ifdef DEBUG
    assert(m_storage.num_elements() > 0);
#endif
  };

  /** Destructor; free memmory */
  ~CoeffMatrix2D() noexcept {
    printf("[DEBUG] CoeffMatrix2D destructor; mem at %p\n", (void*)m_data);
    if (m_data)
      delete[] m_data;
    _capacity = 0;
  }

  /** Copy constructor */
  CoeffMatrix2D(const CoeffMatrix2D &mat) noexcept
      : m_storage(mat.m_storage), m_data(new double[mat.num_elements()]),
        _capacity(mat.num_elements()) {
    std::memcpy(m_data, mat.m_data, sizeof(double) * mat.num_elements());
    printf("[DEBUG] CoeffMatrix2D copy c'tor; mem at: %p size=%ld\n", (void*)m_data, _capacity);
  }

  /** Move constructor */
  CoeffMatrix2D(CoeffMatrix2D &&mat) noexcept
      : m_storage(mat.m_storage), m_data(mat.m_data), _capacity(mat._capacity) {
    mat.m_data = nullptr;
    mat.m_storage.__set_dimensions(0, 0);
    mat._capacity = 0;
    printf("[DEBUG] CoeffMatrix2D move c'tor; mem at: %p size=%ld\n", (void*)m_data, _capacity);
  }

  template <typename T1, typename T2>
  CoeffMatrix2D(_SumProxy<T1, T2> &&sum) noexcept
      : m_storage(sum.rows(), sum.cols()),
        m_data(new double[m_storage.num_elements()]),
        _capacity(m_storage.num_elements()) {
    for (std::size_t i = 0; i < m_storage.num_elements(); i++) {
      m_data[i] = sum.data(i);
    }
    printf("[DEBUG] Called CoeffMatrix2D(_SumProxy<T1, T2> &&sum)\n");
  }
  template <typename T1>
  CoeffMatrix2D(_ScaledProxy<T1> &&fac) noexcept
      : m_storage(fac.rows(), fac.cols()),
        m_data(new double[m_storage.num_elements()]),
        _capacity(m_storage.num_elements()) {
    for (std::size_t i = 0; i < m_storage.num_elements(); i++) {
      m_data[i] = fac.data(i);
    }
    printf("[DEBUG] Called CoeffMatrix2D(_ScaledProxy<T1> &&fac)\n");
  }

  /** (Copy) Assignment operator */
  CoeffMatrix2D &operator=(const CoeffMatrix2D &mat) noexcept {
    if (this != &mat) {
      /* do we need extra capacity ? */
      if (_capacity < mat._capacity) {
        if (m_data)
          delete[] m_data;
        m_data = new double[mat.num_elements()];
        _capacity = mat.num_elements();
      }
      std::memcpy(m_data, mat.m_data, sizeof(double) * mat.num_elements());
      m_storage.__set_dimensions(mat.rows(), mat.cols());
    }
    printf("[DEBUG] Called CoeffMatrix2D &operator=(const CoeffMatrix2D &mat)\n");
    return *this;
  }

  /** Move Assignment operator */
  CoeffMatrix2D &operator=(CoeffMatrix2D &&mat) noexcept {
    m_data = mat.m_data;
    m_storage.__set_dimensions(mat.rows(), mat.cols());
    _capacity = mat._capacity;
    mat.m_data = nullptr;
    mat.m_storage.__set_dimensions(0, 0);
    mat._capacity = 0;
    printf("[DEBUG] Called CoeffMatrix2D &operator=(CoeffMatrix2D &&mat)\n");
    return *this;
  }

  /** Resize (keeping the MatrixStorageType the same) */
  void resize(int rows, int cols) {
    /* do we need to re-allocate ? */
    if (StorageImplementation<S>(rows, cols).num_elements() > _capacity) {
      if (m_data)
        delete[] m_data;
      m_data = new double[StorageImplementation<S>(rows, cols).num_elements()];
      _capacity = StorageImplementation<S>(rows, cols).num_elements();
    } else {
      ;
    }
    m_storage = StorageImplementation<S>(rows, cols);
  }

  /** Sum of two CoeffMatrix2D instances create a proxy instance _SumProxy */
  _SumProxy<CoeffMatrix2D, CoeffMatrix2D>
  operator+(const CoeffMatrix2D &rhs) const noexcept {
    return _SumProxy<CoeffMatrix2D, CoeffMatrix2D>(*this, rhs);
  }

  /** Sum of a CoeffMatrix2D and a Scaled CoeffMatrix2D creates a proxy
   * instance _SumProxy
   */
  template <typename T>
  _SumProxy<CoeffMatrix2D, _ScaledProxy<T>>
  operator+(const _ScaledProxy<T> &rhs) const noexcept {
    return _SumProxy<CoeffMatrix2D, _ScaledProxy<T>>(*this, rhs);
  }

  /** CoeffMatrix2D * Real is a Scaled matrix */
  _ScaledProxy<CoeffMatrix2D> operator*(double factor) const noexcept {
    return _ScaledProxy<CoeffMatrix2D>(*this, factor);
  }

  template <typename T> CoeffMatrix2D &operator+=(const T &rhs) noexcept {
    assert((this->rows() == rhs.rows()) && (this->cols() == rhs.cols()));
    for (std::size_t i = 0; i < m_storage.num_elements(); i++) {
      m_data[i] += rhs.data(i);
    }
    return *this;
  }

}; /* class CoeffMatrix2D */

template <MatrixStorageType S>
inline void swap(CoeffMatrix2D<S> &a, CoeffMatrix2D<S> &b) noexcept {
  a.swap(b);
}

/** Real * CoeffMatrix2D returns a ScaledProxy instance */
template <MatrixStorageType S>
typename CoeffMatrix2D<S>::template _ScaledProxy<CoeffMatrix2D<S>>
operator*(double factor, const CoeffMatrix2D<S> &mat) noexcept {
  return typename CoeffMatrix2D<S>::template _ScaledProxy<CoeffMatrix2D<S>>(
      mat, factor);
}

template <MatrixStorageType S, typename T1, typename T2, typename T3>
typename CoeffMatrix2D<S>::template _SumProxy<
    typename CoeffMatrix2D<S>::template _SumProxy<T1, T2>,
    typename CoeffMatrix2D<S>::template _ScaledProxy<T3>>
operator+(
    const typename CoeffMatrix2D<S>::template _ScaledProxy<T3> &scaled,
    const typename CoeffMatrix2D<S>::template _SumProxy<
        typename CoeffMatrix2D<S>::template _SumProxy<T1, T2>> &sum) noexcept {
  return typename CoeffMatrix2D<S>::template _SumProxy<
      typename CoeffMatrix2D<S>::template _SumProxy<T1, T2>,
      typename CoeffMatrix2D<S>::template _ScaledProxy<T3>>(sum, scaled);
}

} /* namespace dso */

#endif
