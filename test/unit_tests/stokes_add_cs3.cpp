#include "stokes_coefficients.hpp"
#include "iersconst.hpp"
#include <random>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

using namespace dso;


const int MD1 = 180;
const int MO1 = 180;
const int MD2 = 80;
const int MO2 = 55;
constexpr double TOLERANCE = 1e-15;

int main() {

  std::uniform_real_distribution<double> unif(-1e0, 1e0);
  std::default_random_engine re;

  StokesCoeffs cs1(MD1, MO1, iers2010::GMe, iers2010::Re);
  StokesCoeffs cs2(MD2, MO2, iers2010::GMe, iers2010::Re);
  StokesCoeffs res(MD1, MO1, iers2010::GMe, iers2010::Re);

  /* populate cs1 and cs2 with random numbers */
  for (int r=0; r<=MD1; r++) {
    for (int c=0; c<=std::min(r,MO1); c++) {
      cs1.C(r,c) = unif(re);
      cs1.S(r,c) = unif(re);
    }
  }
  for (int r=0; r<=MD2; r++) {
    for (int c=0; c<=std::min(r,MO2); c++) {
      cs2.C(r,c) = unif(re);
      cs2.S(r,c) = unif(re);
    }
  }

  /* res = cs1 + cs2 */
  for (int r=0; r<=MD1; r++) {
    for (int c=0; c<=std::min(r,MO1); c++) {
      if (r>MD2 || c>MO2) {
        res.C(r,c) = cs1.C(r,c);
        res.S(r,c) = cs1.S(r,c);
      } else {
        res.C(r,c) = cs1.C(r,c)+cs2.C(r,c);
        res.S(r,c) = cs1.S(r,c)+cs2.S(r,c);
      }
    }
  }
  
  cs1 += cs2;
  for (int r=0; r<=MD1; r++) {
    for (int c=0; c<=std::min(r,MO1); c++) {
        assert(std::abs(res.C(r,c) - cs1.C(r,c))<TOLERANCE);
        assert(std::abs(res.S(r,c) - cs1.S(r,c))<TOLERANCE);
    }
  }
  
  return 0;
}
