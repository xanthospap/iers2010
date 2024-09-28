#include "coeff_matrix_2d.hpp"
#include <cmath>

using namespace dso;

CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise>
NormalizedLegendre(int n_max, int m_max, double t) {

  assert(n_max > 1 && m_max > 1);
  assert(n_max >= m_max);

  CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> Pnm(n_max + 1,
                                                                 m_max + 1);

  Pnm(0, 0) = 1;
  Pnm(1, 0) = std::sqrt(3) * std::cos(t);
  Pnm(1, 1) = std::sqrt(3) * std::sin(t);

  for (int n = 2; n <= n_max; n++) {
    double l = (2e0 * n + 1e0) / (2e0 * n);
    Pnm(n, n) = std::sin(t) * std::sqrt(l) * Pnm(n - 1, n - 1);
  }

  int max_col = m_max;
  if (n_max == m_max) { max_col--; }
  for (int m = 1; m <= max_col; m++) {
    double l = (2e0 * m + 3e0);
    Pnm(m + 1, m) = std::cos(t) * std::sqrt(l) * Pnm(m, m);
  }

  for (int m = 0; m <= m_max; m++) {
    for (int n = 2; n <= n_max; n++) {
      double l = (2e0 * n - 1e0) * (2e0 * n + 1e0) / (n - m) * (n + m);
      double k = (2e0 * n + 1) * (n + m + 1e0) * (n - m - 1e0) / (n - m) *
                 (n + m) * (2e0 * n - 3e0);
      Pnm(n, m) = std::sqrt(l) * std::cos(t) * Pnm(n - 1e0, m) -
                  std::sqrt(k) * Pnm(n - 2e0, m);
    }
  }
return Pnm;
}
