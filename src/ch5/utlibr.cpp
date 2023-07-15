#include "iers2010.hpp"

namespace {
/* Coefficients of the quasi semidiurnal terms in dUT1, dLOD
 * Source: IERS Conventions (2010), Table 5.1b
 */
constexpr const struct {
  int iarg[6];
  double per, dut1s, dut1c, dlods, dlodc;
} x[] = {{{2, -2, 0, -2, 0, -2}, 0.5377239e0, 0.05e0, -0.03e0, -0.3e0, -0.6e0},
         {{2, 0, 0, -2, -2, -2}, 0.5363232e0, 0.06e0, -0.03e0, -0.4e0, -0.7e0},
         {{2, -1, 0, -2, 0, -2}, 0.5274312e0, 0.35e0, -0.20e0, -2.4e0, -4.1e0},
         {{2, 1, 0, -2, -2, -2}, 0.5260835e0, 0.07e0, -0.04e0, -0.5e0, -0.8e0},
         {{2, 0, 0, -2, 0, -1}, 0.5175645e0, -0.07e0, 0.04e0, 0.5e0, 0.8e0},
         {{2, 0, 0, -2, 0, -2}, 0.5175251e0, 1.75e0, -1.01e0, -12.2e0, -21.3e0},
         {{2, 1, 0, -2, 0, -2}, 0.5079842e0, -0.05e0, 0.03e0, 0.3e0, 0.6e0},
         {{2, 0, -1, -2, 2, -2}, 0.5006854e0, 0.04e0, -0.03e0, -0.3e0, -0.6e0},
         {{2, 0, 0, -2, 2, -2}, 0.5000000e0, 0.76e0, -0.44e0, -5.5e0, -9.6e0},
         {{2, 0, 0, 0, 0, 0}, 0.4986348e0, 0.21e0, -0.12e0, -1.5e0, -2.6e0},
         {{2, 0, 0, 0, 0, -1}, 0.4985982e0, 0.06e0, -0.04e0, -0.4e0, -0.8e0}};
constexpr const int M{sizeof(x) / sizeof(x[0])};
static_assert(M == 11, "Invalid size for quasi semidiurnal terms in utlibr.");
} /* unnamed namespace */

int iers2010::utlibr([[maybe_unused]] const dso::TwoPartDate &mjd,
                     const double *const fargs, double &dut1,
                     double &dlod) noexcept {
  /*
   *         ----------------------------
   *           D E F I N I T I O N S
   *         ----------------------------
   *  iarg   - array defining for each of the 11 trigonometric terms a set
   *           of 6 integer multipliers of the fundamental angular arguments
   *  arg    - vector of the following 6 fundamental arguments used to
   *           compute the angular argument of the trigonometric functions
   *           arg(1:6) = [ GMST+pi, el, elp, f, d, om ]; this vector is
   *           evaluated by the subroutine FUNDARG which is called as an
   *           external subroutine.  Originally evaluated by the subroutine
   *           PMARGS.
   *  period - array of periods of the trigonometric terms of expansion, in
   *           mean solar days; only for a check - not used in computations
   *  dUT1s, dUT1c - sine and cosine coefficients of dUT1, in microseconds
   *  dLODs, dLODc - sine and cosine coefficients of dLOD, in microseconds
   *                 per day
   *  angle  - angular argument of the trigonometric functions
   *           angle = Sum(i=1:6) iarg(i,j)*arg(i), for j=1,11
   */

  /* Compute the harmonic model of dUT1 and dLOD
   * dUT1 and dLOD are set to zero first
   */
  dut1 = dlod = 0e0;

  for (int j = 0; j < M; j++) {
    /* For the jth term of the trigonometric expansion, compute the angular
     * argument angle of sine and cosine functions as a linear integer
     * combination of the 6 fundamental arguments
     */
    double angle = 0e0;
    for (int i = 0; i < 6; i++)
      angle += x[j].iarg[i] * fargs[i];
    angle = dso::anp(angle);
    /* Compute contribution from the jth term of expansion to dUT1 and dLOD */
    const double sa = std::sin(angle);
    const double ca = std::cos(angle);
    dut1 += x[j].dut1s * sa + x[j].dut1c * ca;
    dlod += x[j].dlods * sa + x[j].dlodc * ca;
  }

  /* Finished */
  return 0;
}

int iers2010::utlibr(const dso::TwoPartDate &mjd, double dut1_in,
                     double &dut1_out, double &dlod) noexcept {
  double fargs[6];
  iers2010::utils::fargs(mjd, dut1_in, fargs);
  return utlibr(mjd, fargs, dut1_out, dlod);
}

/*
 * Dear all,
I have a question regarding the UTLIBR.F subroutine. The purpose of the function
is to compute corrections for UT1 (and LOD). However, if we want to be compliant
with the IAU2006 model, the computations performed inside the function (i.e.
angular arguments per tidal constituent) require the value of GMST and hence the
actual UT1 value. The distributed FORTRAN implementation avoids this by not
using the IAU2006 GMST model. What is the recommended way to bypass this
situation? Should one avoid using the IAU2006 GMST model or is it better to use
the  IAU2006 GMST model and use a non-corrected UT1 value (e.g. the ones
published by the IERS)? Thank you very much in advance, Xanthos Papanikolaou
*/
