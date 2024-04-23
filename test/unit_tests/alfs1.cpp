#include "gravity.hpp"
#include <random>
#include <cmath>

int explicit_sh(double lat, dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> &P) {
    /* trigonometric numbers (of third body) */
  const double t = std::sin(lat);
  const double u = std::cos(lat);
  const double u2 = u * u;
  const double t2 = t * t;

  /* compute normalized associated Lagrange polynomials for n=2,3 */
  P(2,0) = std::sqrt(5e0) * 0.5e0 * (3e0 * t2 - 1e0);      // P20
  P(2,1) = std::sqrt(5.e0 / 3.e0) * 3.e0 * t * u;          // P21
  P(2,2) = std::sqrt(5.e0 / 12.e0) * 3.e0 * u2;            // P22
  P(3,0) = (.5e0 * t * std::sqrt(7e0)) * (5e0 * t2 - 3e0); // P30
  P(3,1) =
    (3e0 / 2e0) * (5e0 * t2 - 1e0) * u * std::sqrt((7e0) / 6e0); // P31
  P(3,2) = 15e0 * std::sqrt(7e0 / 60e0) * t * u2;      // P32
  P(3,3) = 15e0 * u2 * u * std::sqrt(14e0 / 720e0); // P33

  return 0;
}

void dprint(dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> &P) {
  printf("%+15.9f %+15.9f %+15.9f\n", P(2,0), P(2,1), P(2,2));
  printf("%+15.9f %+15.9f %+15.9f %+15.9f\n", P(3,0), P(3,1), P(3,2), P(3,3));
}

constexpr const int num_tests=100;
constexpr const double PRECISION = 1e-14;
int main() {
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> P1(5, 5),
      P2(5, 5);
   std::uniform_real_distribution<double> unif(-M_PI, M_PI);
   std::default_random_engine re;

  for (int i=0; i<num_tests; i++) {
    const double theta = unif(re);
    explicit_sh(theta, P1);
    // dprint(P1);

    dso::normalised_associated_legendre_functions(theta, 4, 4, P2);
    // dprint(P2);

    assert(std::abs(P1(2,0) - P2(2,0))<PRECISION);
    assert(std::abs(P1(2,1) - P2(2,1))<PRECISION);
    assert(std::abs(P1(2,2) - P2(2,2))<PRECISION);
    assert(std::abs(P1(3,0) - P2(3,0))<PRECISION);
    assert(std::abs(P1(3,1) - P2(3,1))<PRECISION);
    assert(std::abs(P1(3,2) - P2(3,2))<PRECISION);
    assert(std::abs(P1(3,3) - P2(3,3))<PRECISION);
    
    dso::normalised_associated_legendre_functions(theta, 3, 3, P2);
    // dprint(P2);

    assert(std::abs(P1(2,0) - P2(2,0))<PRECISION);
    assert(std::abs(P1(2,1) - P2(2,1))<PRECISION);
    assert(std::abs(P1(2,2) - P2(2,2))<PRECISION);
    assert(std::abs(P1(3,0) - P2(3,0))<PRECISION);
    assert(std::abs(P1(3,1) - P2(3,1))<PRECISION);
    assert(std::abs(P1(3,2) - P2(3,2))<PRECISION);
    assert(std::abs(P1(3,3) - P2(3,3))<PRECISION);
  }

  return 0;
}
