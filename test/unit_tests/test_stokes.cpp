#include "stokes_coefficients.hpp"
#include <cassert>
#include <cstdio>
#include <random>
#include <typeinfo>

using namespace dso;
constexpr const int N = 180;
constexpr const int M = 120;
std::uniform_real_distribution<double> unif(-1e0, 1e0);
std::default_random_engine re;

void fillRandom(StokesCoeffs &cs) {
  for (int m = 0; m <= cs.max_order(); m++) {
    for (int n = m; n <= cs.max_degree(); n++) {
      cs.C(n, m) = unif(re);
      cs.S(n, m) = unif(re);
    }
  }
  return;
}

int main() {
  constexpr const int numtests = 10;
  StokesCoeffs a2(N, M, 0e0, 0e0), a3(N, M, 0e0, 0e0), a4(N, M, 0e0, 0e0);

  assert(a2.normalized() == a3.normalized());
  assert(a3.normalized() == a4.normalized());

  double dummy = 0e0;

  for (int i = 0; i < numtests; i++) {
    fillRandom(a2);
    fillRandom(a3);
    fillRandom(a4);
    /* Note 1 */
    auto a1 = 0.1 * a2 + 0.2 * a3 + (-.1) * a4;
    dummy += a1.C(0, 0) + a1.S(1, 1);
    printf("\t%.6f <- %.9f * %.3f + %.9f * %.3f + %.9f * %.3f\n", a1.C(3, 3),
           a2.C(3, 3), .1, a3.C(3, 3), .2, a4.C(3, 3), -.1);
    printf("Size of a1 = (%dx%d), type is %s\n", a1.max_degree(),
           a1.max_order(), typeid(a1).name());
  }
  printf("Exit loop 1\n");

/* Note 1:
 * What kind of instance would you say a1 is ?
 * You are probably wrong; it is not a StokesCoeffs instance, it is an 
 * instance of type StokesCoeffs::_SumProxy<ScaledProxy...>
 *
 * In the following loop however, a1 is indeed a StokesCoeffs instance.
 */
  for (int i = 0; i < numtests; i++) {
    fillRandom(a2);
    fillRandom(a3);
    fillRandom(a4);
    /* Note 1 */
    StokesCoeffs a1 = 0.1 * a2 + 0.2 * a3 + (-.1) * a4;
    dummy += a1.C(0, 0) + a1.S(1, 1);
    printf("\t%.6f <- %.9f * %.3f + %.9f * %.3f + %.9f * %.3f\n", a1.C(3, 3),
           a2.C(3, 3), .1, a3.C(3, 3), .2, a4.C(3, 3), -.1);
    printf("Size of a1 = (%dx%d), type is %s\n", a1.max_degree(),
           a1.max_order(), typeid(a1).name());
  }
  printf("Exit loop 1.1\n");
  
  { /* should NOT produce an error */
    StokesCoeffs er(N-1,M-1,0e0,0e00);
    a2 += er;
    dummy += a2.C(0, 0) + a2.S(1, 1);
  }
  printf("Exit loop 2\n");
  { /* should NOT produce an error */
    StokesCoeffs b2(N - 1, M - 1, 0e0, 0e0), b3(N - 1, M - 1, 0e0, 0e0),
        b4(N - 1, M - 1, 0e0, 0e0);
    StokesCoeffs b(N, M, 0e0, 0e0);
    b += 0.1 * b2 + 0.2 * b3 + (-.1) * b4;
    dummy += b.C(0, 0) + b.S(1, 1);
  }
  printf("Exit loop 3\n");

  //{ /* should produce an error */
  //  StokesCoeffs er(N,M,1e0,0e00);
  //  auto a1 = 0.1 * er + 0.2 * a2;
  //  dummy += a1.C(0, 0) + a1.S(1, 1);
  //}
  //{ /* should produce an error */
  //  StokesCoeffs er(N,M,0e0,1e00);
  //  auto a1 = 0.1 * er + 0.2 * a2;
  //  dummy += a1.C(0, 0) + a1.S(1, 1);
  //}
  //{ /* should produce an error */
  //  StokesCoeffs er(N-1,M,0e0,0e00);
  //  auto a1 = 0.1 * er + 0.2 * a2;
  //  dummy += a1.C(0, 0) + a1.S(1, 1);
  //}

  printf("Dummy: %.5f\n", dummy);
  return 0;
}
