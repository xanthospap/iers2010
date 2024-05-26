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

  double f[5];
  for (int i=0;i<5;i++) f[i] = unif(re);

  cs.clear();
  for (int j=0; j<5; j++) {
    cs.Cnm() += f[j] * cs1.Cnm() + (f[j]*.1) * cs2.Cnm();
    cs.Snm() += f[j] * cs1.Snm() + (f[j]*.1) * cs2.Snm();
  }

  /* validation */
  for (int c=0; c<=180; c++) {
    for (int r=c; r<=c; r++) {
      double cnm = 0e0;
      double snm = 0e0;
      for (int j=0;j<5;j++) {
        cnm += f[j] * cs1.C(r,c) + (f[j]*.1) * cs2.C(r,c);
        snm += f[j] * cs1.S(r,c) + (f[j]*.1) * cs2.S(r,c);
      }
      assert( cs.C(r,c) == cnm );
      assert( cs.S(r,c) == snm );
    }
  }

  return 0;
}
