#include "stokes_coefficients.hpp"
#include "iersconst.hpp"
#include <random>
#include <cassert>

using namespace dso;


const int MAX_DEGREE = 180;
const int MAX_ORDER = 180;
int main() {

  std::uniform_real_distribution<double> unif(-1e0, 1e0);
  std::default_random_engine re;

  StokesCoeffs cs (MAX_DEGREE, MAX_ORDER, iers2010::GMe, iers2010::Re);
  StokesCoeffs cs1(MAX_DEGREE, MAX_ORDER, iers2010::GMe, iers2010::Re);
  StokesCoeffs cs2(MAX_DEGREE, MAX_ORDER, iers2010::GMe, iers2010::Re);

  /* populate cs1 and cs2 with random numbers */
  for (int c=0; c<=180; c++) {
    for (int r=c; r<=c; r++) {
      cs1.C(r,c) = unif(re);
      cs1.S(r,c) = unif(re);
      cs2.C(r,c) = unif(re);
      cs2.S(r,c) = unif(re);
    }
  }

  cs.clear();
  cs.Cnm() = 1.5e0 * cs1.Cnm() + 2.5e0 * cs2.Cnm();
  cs.Snm() = 1.5e0 * cs1.Snm() + 2.5e0 * cs2.Snm();

  /* validation */
  for (int c=0; c<=180; c++) {
    for (int r=c; r<=c; r++) {
      assert( cs.C(r,c) == cs1.C(r,c) * 1.5e0 + cs2.C(r,c) * 2.5e0 );
      assert( cs.S(r,c) == cs1.S(r,c) * 1.5e0 + cs2.S(r,c) * 2.5e0 );
    }
  }

  return 0;
}
