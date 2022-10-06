#ifdef USE_EIGEN
#include <random>
#include <cassert>
#include <limits>
#include <cmath>
#include "iau.hpp"

bool equal(const dso::Mat3x3& r, const Eigen::Matrix<double, 3, 3> &re,
           double tolerance = std::numeric_limits<double>::epsilon()) {
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      if (std::abs(re(i,j) - r(i,j)) > tolerance) return false;
    }
  }
  return true;
}

template<typename T>
void printm(const T& mat, const char *title=nullptr) {
  if (title) printf("%s\n", title);
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      printf(" %+.6f ", mat(i,j));
    }
    printf("\n");
  }
  return;
}

std::uniform_real_distribution<double> unif(-2e0*M_PI, 2e0*M_PI);
std::default_random_engine re;
int main() {

  constexpr const int num_tests = 10;
  
  printf("Testing function pom00 ...\n");
  for (int t=0; t<num_tests; t++) {
    const double x = unif(re)  * 1e2;
    const double y = unif(re)  * 1e2;
    const double a = unif(re)  * 1e2;
    const auto m1 = iers2010::sofa::pom00_e(x, y, a);
    const auto m2 = iers2010::sofa::pom00(x, y, a);
    if (!equal(m2,m1,1e-15)) {
      printf("Test nr: %d\n", t);
      printm(m1, "Eigen ->");
      printm(m2, "Sofa  ->");
    }
    assert(equal(m2,m1,1e-15));
  }
  
  printf("Testing function c2ixys...\n");
  for (int t=0; t<num_tests; t++) {
    const double x = unif(re)  * 1e2;
    const double y = unif(re)  * 1e2;
    const double a = unif(re)  * 1e2;
    const auto m1 = iers2010::sofa::c2ixys_e(x, y, a);
    const auto m2 = iers2010::sofa::c2ixys(x, y, a);
    if (!equal(m2,m1,1e-15)) {
      printf("Test nr: %d\n", t);
      printm(m1, "Eigen ->");
      printm(m2, "Sofa  ->");
    }
    assert(equal(m2,m1,1e-15));
  }

  printf("All tests ok!\n");
}

#endif
