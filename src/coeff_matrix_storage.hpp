/** @file
 * Declerations of MatrixStorageType, i.e. indexing methods dependent on
 * the way a 2D matrix is stored in memory. We are not interested here in
 * the actual memmory (allocating/de-allocationg etc), but only on indexing.
 */

#ifndef __COMPACT_2D_SIMPLE_MATRIX_STORAGE_HPP__
#define __COMPACT_2D_SIMPLE_MATRIX_STORAGE_HPP__

#include <algorithm>
#include <cstddef>
#include <cstdio>
#ifdef DEBUG
#include <cassert>
#include <cstdio>
#endif

namespace dso {

/** Enum class to describe storage type of a 2-d matrix */
enum class MatrixStorageType : char {
  RowWise,             /** Row-Wise storage */
  ColumnWise,          /** Column-Wise storage */
  Trapezoid,           /** A trapezoid matrix with row-wise storage */
  LwTriangularRowWise, /** Lower triangular, stored Row-Wise */
  LwTriangularColWise  /** Lower triangular, stored Col-Wise */
};                     /* MatrixStorageType */

/** @brief implementation details depending on storage type, aka
 *    MatrixStorageType
 */
template <MatrixStorageType S> class StorageImplementation {};

/** @brief Implementation details for a 2-d lower triangular matrix, holding
 *        data in a Row-Wise fashion.
 *
 * Here we are not interested on the actual data of the matrix, but only the
 * indexing implementation of the matrix. These are always considered as
 * rectangular matrices, i.e. the number of rows is the same as the number of
 * columns.
 */
template <>
class StorageImplementation<MatrixStorageType::LwTriangularRowWise> {
private:
  /** num of rows (= number of columns) */
  int rows;

  /** number of elements in a Lower-Triangular square matrix with size=d */
  constexpr std::size_t _size(int d) const noexcept { return d * (d + 1) / 2; }

public:
  /** Basic stride/dimension */
  static constexpr const bool isRowMajor = true;
  static constexpr const bool isColMajor = false;

  /** Constructor; not interested in number of cols */
  constexpr StorageImplementation(int r, [[maybe_unused]] int c) noexcept
      : rows(r)
  {
#ifdef DEBUG
    assert(r==c);
#endif
  };

  /** (Re-) set dimensions */
  void __set_dimensions(int _rows, [[maybe_unused]] int _cols) noexcept {
    rows = _rows;
#ifdef DEBUG
    assert(_rows==_cols);
#endif
  }

  /** @brief Compute number of elements stored */
  constexpr std::size_t num_elements() const noexcept { return _size(rows); }

  /** get number of rows */
  constexpr int nrows() const noexcept { return rows; }

  /** get number of cols */
  constexpr int ncols() const noexcept { return rows; }

  /** @brief Index of (beggining of) row
   *
   * Return the offset from the begining of the data array, given a row
   * number. First row is row 0 (NOT row 1).
   * That means that if the data is stored in an array e.g.
   *   double *data = new double[num_pts];
   *   double *row_3 = data[0] + slice(2);
   * will point to the first (0) element of the third row.
   */
  constexpr int slice(int row) const noexcept { return _size(row); }
  
  /** @brief Index of (beggining of) row and number of elements in row.
   *
   * Return the offset from the begining of the data array, given a row
   * number. First row is row 0 (NOT row 1).
   * The input parameter \p num_elements will be set to the nu,ber of elements 
   * stored in the row requested.
   */
  constexpr int slice(int row, int &num_elements) const noexcept {
    num_elements = row + 1;
    return _size(row); 
  }

  /** @brief Number of slices, i.e. number of rows */
  constexpr int num_slices() const noexcept { return rows; }

  /** @brief Index of element (row, column) in the data array.
   *  E.g. data[element_offset(1,2)] will return the element in the second
   *  row, and third column.
   */
  constexpr int index(int row, int column) const noexcept {
    return slice(row) + column;
  }
}; /* StorageImplementation<MatrixStorageType::LwTriangularRowWise> */

/** @brief Implementation details for a 2-d lower triangular matrix, holding
 *        data in a Col-Wise fashion.
 *
 * Here we are not interested on the actual data of the matrix, but only the
 * indexing implementation of the matrix. These are always considered as
 * rectangular matrices, i.e. the number of rows is the same as the number of
 * columns.
 */
template <>
class StorageImplementation<MatrixStorageType::LwTriangularColWise> {
private:
  /** num of rows (= number of columns) */
  int rows;

public:
  /** Basic stride/dimension */
  static constexpr const bool isRowMajor = false;
  static constexpr const bool isColMajor = true;

  /** Constructor; not interested in number of cols */
  constexpr StorageImplementation(int r, [[maybe_unused]] int c) noexcept
      : rows(r)
      {
#ifdef DEBUG
    assert(r==c);
#endif
      }

  /** @brief Compute number of elements stored */
  constexpr std::size_t num_elements() const noexcept {
    return rows * (rows + 1) / 2;
  }

  /** (Re-) set dimensions */
  void __set_dimensions(int _rows, [[maybe_unused]] int _cols) noexcept {
    rows = _rows;
#ifdef DEBUG
    assert(_rows==_cols);
#endif
  }

  /** number of rows */
  constexpr int nrows() const noexcept { return rows; }

  /** number of columns */
  constexpr int ncols() const noexcept { return rows; }

  /** @brief Index/offset of given column.
   *
   * Return the offset from the begining of the data array, given a column
   * number.
   * First column is column 0 (NOT row 1).
   * That means that if the data is stored in an array e.g.
   *   double *data = new double[num_pts];
   *   double *col_3 = data[0] + slice(2);
   * will point to the first (0) element of the third column.
   */
  constexpr int slice(int col) const noexcept {
    return col * rows - col * (col - 1) / 2;
  }
  
  constexpr int slice(int col, int &num_elements) const noexcept {
    num_elements = col + 1;
    return slice(col);
  }
  
  /** @brief Number of slices, i.e. number of cols */
  constexpr int num_slices() const noexcept { return rows; }

  /** @brief Index of element (row, column) in the data array.
   *  E.g. data[element_offset(1,2)] will return the element in the second
   *  row, and third column.
   */
  constexpr int index(int row, int col) const noexcept {
    return slice(col) + (row - col);
  }
}; /* StorageImplementation<MatrixStorageType::LwTriangularColWise> */

/** @brief Implementation details for a 2-d trapezoid matrix, holding data in
 *        a Row-Wise fashion.
 *
 * In case rows == columns, this is actually a triangular matrix.
 * Here we are not interested on the actual data of the matrix, but only the
 * indexing implementation of the matrix.
 */
template <> class StorageImplementation<MatrixStorageType::Trapezoid> {
private:
  /** number of rows */
  int rows;
  /** number of columns */
  int cols;

public:
  /** Constructor given number of rows and num of columns */
  constexpr StorageImplementation(int r, int c) noexcept : rows(r), cols(c){};

  /** number of rows */
  constexpr int nrows() const noexcept { return rows; }

  /** number of columns */
  constexpr int ncols() const noexcept { return cols; }

  /** (Re-) set dimensions */
  void __set_dimensions(int _rows, int _cols) noexcept {
    rows = _rows;
    cols = _cols;
  }

  /** @brief Compute number of elements stored */
  constexpr std::size_t num_elements() const noexcept {
    if (rows == cols) {
      std::size_t n = rows;
      return (n * (n + 1)) / 2;
    } else {
      std::size_t n = 0;
      for (int r = 0; r < rows; r++)
        n += std::min(r + 1, cols);
      return n;
    }
  }

  /** @brief Number of elements (data points) stored for a given row */
  constexpr int pts_in_row(int row) const noexcept {
    return std::min(row + 1, cols);
  }

  /** @brief Index of (beggining of) row.
   *
   * Return the offset from the begining of the data array, given a row
   * number. First row is row 0 (NOT row 1).
   * That means that if the data is stored in an array e.g.
   *   double *data = new double[num_pts];
   *   double *row_3 = data[0] + slice(2);
   * will point to the first (0) element of the third row.
   */
  constexpr int slice(int row) const noexcept {
    int offset = 0;
    for (int i = 0; i < row; i++)
      offset += pts_in_row(i);
    return offset;
  }
  
  constexpr int slice(int row, int &num_elements) const noexcept {
    num_elements = pts_in_row(row);
    return slice(row);
  }
  
  /** @brief Number of slices, i.e. number of rows */
  constexpr int num_slices() const noexcept { return rows; }

  /** @brief Index of element (row, column) in the data array.
   * E.g. data[element_offset(1,2)] will return the element in the second row,
   *  and third column.
   */
  constexpr int index(int row, int column) const noexcept {
    return (column <= rows) ? (slice(row) + column) : (slice(column) + row);
  }

}; /* StorageImplementation<MatrixStorageType::Trapezoid> */

/** @brief Implementation details for a 2-d dense matrix, holding data in
 *        a Row-Wise fashion.
 *  Here we are not interested on the actual data of the matrix, but only the
 *  indexing implementation of the matrix.
 */
template <> struct StorageImplementation<MatrixStorageType::RowWise> {
private:
  /** number of rows */
  int rows;
  /** number of columns */
  int cols;

public:
  /** Basic stride/dimension */
  static constexpr const bool isRowMajor = true;
  static constexpr const bool isColMajor = false;

  /** Constructor given number of rows and number of columns */
  constexpr StorageImplementation(int r, int c) noexcept : rows(r), cols(c){};

  /** get number of rows */
  constexpr int nrows() const noexcept { return rows; }

  /* get number of columns */
  constexpr int ncols() const noexcept { return cols; }

  /** (Re-) set dimensions */
  void __set_dimensions(int _rows, int _cols) noexcept {
    rows = _rows;
    cols = _cols;
  }

  /** @brief Number of elements in matrix */
  constexpr std::size_t num_elements() const noexcept { return rows * cols; }

  /** @brief Index of element (row, column) in the data array.
   *  E.g. data[element_offset(1,2)] will return the element in the second
   *  row, and third column.
   */
  constexpr int index(int row, int column) const noexcept {
    return row * cols + column;
  }

  /** @brief Index/offset of given row.
   *
   * Return the offset from the begining of the data array, given a row
   * number.
   * First row is row 0 (NOT row 1).
   * That means that if the data is stored in an array e.g.
   *   double *data = new double[num_pts];
   *   double *row_3 = data[0] + slice(2);
   * will point to the first (0) element of the third row.
   */
  constexpr int slice(int row) const noexcept { return row * cols; }
  
  constexpr int slice(int row, int &num_elements) const noexcept {
    num_elements = row + 1;
    return slice(row); 
  }
  
  /** @brief Number of slices, i.e. number of rows */
  constexpr int num_slices() const noexcept { return rows; }
}; /* StorageImplementation<MatrixStorageType::RowWise> */

/** @brief Implementation details for a 2-d dense matrix, holding data in
 *        a Column-Wise fashion.
 * Here we are not interested on the actual data of the matrix, but only the
 * indexing implementation of the matrix.
 */
template <> class StorageImplementation<MatrixStorageType::ColumnWise> {
private:
  /** number of rows */
  int rows;
  /** number of columns */
  int cols;

public:
  /** Basic stride/dimension */
  static constexpr const bool isRowMajor = false;
  static constexpr const bool isColMajor = true;

  /** Constructor given number of rows and number of columns */
  constexpr StorageImplementation(int r, int c) noexcept : rows(r), cols(c){};

  /** get number of rows */
  constexpr int nrows() const noexcept { return rows; }

  /** get number of columns */
  constexpr int ncols() const noexcept { return cols; }

  /** (Re-) set dimensions */
  void __set_dimensions(int _rows, int _cols) noexcept {
    rows = _rows;
    cols = _cols;
  }

  /** @brief Number of elements in matrix */
  constexpr std::size_t num_elements() const noexcept { return rows * cols; }

  /** @brief Index of element (row, column) in the data array.
   *  E.g. data[element_offset(1,2)] will return the element in the second
   *  row, and third column.
   */
  constexpr int index(int row, int column) const noexcept {
    return column * rows + row;
  }

  /** @brief Index/offset of given column.
   *
   * Return the offset from the begining of the data array, given a column
   * number.
   * First column is column 0 (NOT row 1).
   * That means that if the data is stored in an array e.g.
   *   double *data = new double[num_pts];
   *   double *col_3 = data[0] + slice(2);
   * will point to the first (0) element of the third column.
   */
  constexpr int slice(int col) const noexcept { return col * rows; }
  
  constexpr int slice(int col, int &num_elements) const noexcept {
    num_elements = ncols();
    return slice(col); 
  }

  constexpr int num_slices() const noexcept { return ncols(); }

}; /* StorageImplementation<MatrixStorageType::ColumnWise> */

} /* namespace dso */

#endif
