/** @file
 * Naive implementation of 2-Dimensional matrices; note that this
 * implementation targets 2-dimensional matrices that are only meant to store
 * data NOT perform arithmetic operations.
 */

#ifndef __COMPACT_2D_SIMPLE_MATRIX_HPP__
#define __COMPACT_2D_SIMPLE_MATRIX_HPP__

#include "coeff_matrix_storage.hpp"
#include <cstring>
#ifdef DEBUG
#include <cassert>
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

public:
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

#ifdef DEBUG
  double *data() noexcept { return m_data; }
#endif

  /** Constructor using number of rows and columns; for some
   * MatrixStorageType's, the number of columns may not be needed.
   */
  CoeffMatrix2D(int rows, int cols) noexcept
      : m_storage(rows, cols), m_data(new double[m_storage.num_elements()]) {
#ifdef DEBUG
    assert(m_storage.num_elements() > 0);
#endif
  };

  /** Destructor; free memmory */
  ~CoeffMatrix2D() noexcept {
    if (m_data)
      delete[] m_data;
  }

  /** Copy constructor */
  CoeffMatrix2D(const CoeffMatrix2D &mat) noexcept
      : m_storage(mat.m_storage), m_data(new double[mat.num_elements()]) {
    std::memcpy(m_data, mat.m_data, sizeof(double) * mat.num_elements());
  }

  /** Move constructor */
  CoeffMatrix2D(CoeffMatrix2D &&mat) noexcept
      : m_storage(mat.m_storage), m_data(mat.m_data) {
    mat.m_data = nullptr;
    mat.m_storage.__set_dimensions(0, 0);
  }

  /** (Copy) Assignment operator */
  CoeffMatrix2D &operator=(const CoeffMatrix2D &mat) noexcept {
    if (this != &mat) {
      if (m_data && (m_storage.num_elements() != mat.num_elements())) {
        delete[] m_data;
        m_data = new double[mat.num_elements()];
      }
      std::memcpy(m_data, mat.m_data, sizeof(double) * mat.num_elements());
      m_storage.__set_dimensions(mat.rows(), mat.cols());
    }
    return *this;
  }

  /** Move Assignment operator */
  CoeffMatrix2D &operator=(CoeffMatrix2D &&mat) noexcept {
    m_data = mat.m_data;
    m_storage.__set_dimensions(mat.rows(), mat.cols());
    mat.m_data = nullptr;
    mat.m_storage.__set_dimensions(0, 0);
    return *this;
  }

  /** Resize (keeping the MatrixStorageType the same) */
  void resize(int rows, int cols) {
    /* de we need to re-allocate ? */
    if (m_storage.num_elements() !=
        StorageImplementation<S>(rows, cols).num_elements()) {
      if (m_data)
        delete[] m_data;
      m_data = new double[StorageImplementation<S>(rows, cols).num_elements()];
    }
    m_storage = StorageImplementation<S>(rows, cols);
  }

}; /* class CoeffMatrix2D */

} /* namespace dso */

#endif
