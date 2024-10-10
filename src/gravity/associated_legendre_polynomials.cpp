#include "gravity.hpp"
#include <cmath>
#include "../coeff_matrix_2d.hpp"

// A function that calculates the normalized associated Legendre Polynomials
// using a recursive algorithm
dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise>
dso::normalised_alfs(int n_max, int m_max, double t) {

  assert(n_max > 1 && m_max > 1);
  assert(n_max >= m_max);

  // Initialize a 2-d lower triangular matrix to store the values of the
  // Legendre Polynomials
  CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> Pnm(n_max + 1,
                                                                 m_max + 1);


  double sint = std::sin(t);
  double cost = std::cos(t);
  
  // Initialize the base values in order to start the algorithm
  Pnm(0, 0) = 1;
  Pnm(1, 0) = std::sqrt(3e0) * cost;
  Pnm(1, 1) = std::sqrt(3e0) * sint;
  
  // The Sectorial (diagonal) values and the off - diagonal values are
  // calculated first and used as seed values for the recursion 
  // Sectorial (i.e. n = m) Pn,n(t)
  for (int n = 2; n <= n_max; n++) {
    double l = (2e0 * n + 1e0) / (2e0 * n);
    Pnm(n, n) = sint * std::sqrt(l) * Pnm(n - 1, n - 1);
  }

  // Off - Diagonal (Degree = Order + 1)
  int max_col = m_max;
  if (n_max == m_max) {
    max_col--;
  }
  for (int m = 1; m <= max_col; m++) {
    double l = (2e0 * m + 3e0);
    Pnm(m + 1, m) = cost * std::sqrt(l) * Pnm(m, m);
  }

  // Non - Sectorial (i.e. n > m) Pn,m(t) from previously computed values
  for (int m = 0; m <= m_max; m++) {
    for (int n = m + 2; n <= n_max; n++) {
      double a = (2e0 * n - 1e0) * (2e0 * n + 1e0) / ((n - m) * (n + m));
      double b = (2e0 * n + 1e0) * (n + m - 1e0) * (n - m - 1e0) / ((n - m) *
                 (n + m) * (2e0 * n - 3e0));
      Pnm(n, m) = std::sqrt(a) * cost * Pnm(n - 1, m) -
                  std::sqrt(b) * Pnm(n - 2, m);
    }
  }
  return Pnm;
}
