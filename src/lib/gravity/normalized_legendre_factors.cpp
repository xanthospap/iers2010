#include "gravity.hpp"
#include <array>
#include <cmath>

/* Max size for ALF factors; if degree is more than this (-2), then it must
 * be augmented. For now, OK
 */
constexpr const int MAX_SIZE_FOR_ALF_FACTORS = 201;

/** @brief  Precompute and store recurrence coefficients (f1 and f2) for
 * efficiently evaluating associated Legendre functions.
 *
 * Additionally the function will precompute and store (array) coefficients f3,
 * which are needed to compute acceleration.
 *
 * Max (n,m) = MAX_SIZE_FOR_ALF_FACTORS, i.e
 * n = [0, ... MAX_SIZE_FOR_ALF_FACTORS-1]
 * m = [0, ... MAX_SIZE_FOR_ALF_FACTORS-1]
 *
 * These factors are used in recursive formulas that generate the Legendre
 * functions (and therefore the spherical harmonic basis functions) without
 * recomputing expensive square roots repeatedly.
 */
struct NormalizedLegendreFactors {
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> f1;
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> f2;
  std::array<double, MAX_SIZE_FOR_ALF_FACTORS> f3;

  NormalizedLegendreFactors() noexcept
      : f1(MAX_SIZE_FOR_ALF_FACTORS, MAX_SIZE_FOR_ALF_FACTORS),
        f2(MAX_SIZE_FOR_ALF_FACTORS, MAX_SIZE_FOR_ALF_FACTORS) {
    constexpr const int N = MAX_SIZE_FOR_ALF_FACTORS;
    f1.fill_with(0e0);
    f2.fill_with(0e0);

    /* factors for the recursion ((2n+1)/2n)^(1/2) */
    f1(1, 1) = std::sqrt(3e0);
    for (int n = 2; n < N; n++) {
      f1(n, n) = std::sqrt((2e0 * n + 1e0) / (2e0 * n));
    }

    /* factors for the recursion */
    for (int m = 0; m < N - 1; m++) {
      for (int n = m + 1; n < N; n++) {
        const double f =
            (2e0 * n + 1e0) / static_cast<double>((n + m) * (n - m));
        /* f1_nm = B_nm */
        f1(n, m) = std::sqrt(f * (2e0 * n - 1e0));
        /* f2_nm = B_nm / Bn-1m */
        f2(n, m) =
            -std::sqrt(f * (n - m - 1e0) * (n + m - 1e0) / (2e0 * n - 3e0));
      }
    }

    /* factors for acceleration */
    for (int n = 0; n < N; n++)
      f3[n] = std::sqrt((double)(2 * n + 1) / (2 * n + 3));
  }
}; /* NormalizedLegendreFactors */