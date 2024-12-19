#include "fundarg.hpp"
#include "iau.hpp"
#include <array>
#include <cmath>
#include <cstdint>

namespace {
/* A struct to hold terms for the expansion of the quantity s(t)+XY/2,
 * using the IAU 2006 precession and IAU 2000A_R06 nutation.
 * See IERS2010, Section 5.5.6 and table 5.2c
 */
struct CioSeriesData {
  /* C_{s,j})_i */
  double csj;
  /* C_{c,j})_i */
  double ccj;
  /* l    l'   F    D   Om L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A */
  int8_t f[14];
  int has_planetary_args() const noexcept {
    for (int i=5; i<14; i++) if (f[i]) return 1;
    return 0;
  }
  /* Find argument given the full list of luni-solar and planetary arguments.
   * That is: ARGUMENT = Σ(i=0,14) N_i * F_i
   * where F_i are given (in the fargs parameter) and N_i are stored in the
   * f[14] member variable.
   */
  double arg(const double *const fargs) const noexcept {
    double arg = 0e0;
    arg = fargs[0] * f[0] + fargs[1] * f[1] + fargs[2] * f[2] +
          fargs[3] * f[3] + fargs[4] * f[4];
    if (has_planetary_args()) {
      for (int i = 5; i < 14; i++)
        arg += fargs[i] * f[i];
    }
    return arg;
  }
}; /* CioSeriesData */

/* Expansion for s(t), series; j=0, Number of terms = 38
 * Extracted from Table 5.2c
 * Note that we are storing elements in reverse order (i.e. from i=38 to 1),
 * to traverse the smallest terms first.
 */
constexpr const std::array<CioSeriesData, 38> C0Series = {
    {{/*   33*/ -1.100000000000e-01,
      +0.000000000000e+00,
      {+1, +0, -2, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   32*/ -1.100000000000e-01,
      +0.000000000000e+00,
      {+1, +0, -2, +0, -3, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   31*/ +1.100000000000e-01,
      +0.000000000000e+00,
      {+0, +0, +2, -2, +4, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   30*/ -1.300000000000e-01,
      +0.000000000000e+00,
      {+0, +0, +4, -2, +4, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   29*/ -1.400000000000e-01,
      +0.000000000000e+00,
      {+1, +0, +0, -2, -1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   28*/ -1.400000000000e-01,
      +0.000000000000e+00,
      {+1, +0, +0, -2, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   27*/ +1.400000000000e-01,
      +0.000000000000e+00,
      {+0, +1, +2, -2, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   26*/ +1.400000000000e-01,
      +0.000000000000e+00,
      {+2, +0, -2, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   25*/ -1.500000000000e-01,
      +0.000000000000e+00,
      {+0, +0, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   24*/ +1.000000000000e-01,
      -5.000000000000e-02,
      {+0, +0, +0, +0, +0, +0, +8, -13, +0, +0, +0, +0, +0, -1}},
     {/*   23*/ -1.800000000000e-01,
      +0.000000000000e+00,
      {+0, +1, -2, +2, -1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   22*/ -1.900000000000e-01,
      +0.000000000000e+00,
      {+0, +1, -2, +2, -3, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   21*/ +2.100000000000e-01,
      +0.000000000000e+00,
      {+0, +0, +2, -2, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   20*/ -2.600000000000e-01,
      +0.000000000000e+00,
      {+1, +0, +2, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   19*/ -2.700000000000e-01,
      +0.000000000000e+00,
      {+1, +0, +2, +0, +3, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   18*/ -2.800000000000e-01,
      +0.000000000000e+00,
      {+0, +0, +2, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   17*/ -3.200000000000e-01,
      +0.000000000000e+00,
      {+0, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   16*/ +2.400000000000e-01,
      +1.200000000000e-01,
      {+0, +0, +1, -1, +1, +0, -8, +12, +0, +0, +0, +0, +0, +0}},
     {/*   15*/ -3.600000000000e-01,
      +0.000000000000e+00,
      {+0, +0, +4, -4, +4, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   14*/ -4.500000000000e-01,
      +0.000000000000e+00,
      {+0, +1, +2, -2, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   13*/ -4.600000000000e-01,
      +0.000000000000e+00,
      {+0, +1, +2, -2, +3, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   12*/ +6.300000000000e-01,
      +0.000000000000e+00,
      {+1, +0, +0, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   11*/ +6.300000000000e-01,
      +0.000000000000e+00,
      {+1, +0, +0, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   10*/ +1.260000000000e+00,
      +1.000000000000e-02,
      {+0, +1, +0, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*    9*/ +1.410000000000e+00,
      +1.000000000000e-02,
      {+0, +1, +0, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*    8*/ +1.720000000000e+00,
      +0.000000000000e+00,
      {+0, +0, +0, +0, +3, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*    7*/ -1.980000000000e+00,
      +0.000000000000e+00,
      {+0, +0, +2, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*    6*/ -2.020000000000e+00,
      +0.000000000000e+00,
      {+0, +0, +2, +0, +3, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*    5*/ +4.570000000000e+00,
      +0.000000000000e+00,
      {+0, +0, +2, -2, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*    4*/ -1.121000000000e+01,
      -1.000000000000e-02,
      {+0, +0, +2, -2, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*    3*/ -1.175000000000e+01,
      -1.000000000000e-02,
      {+0, +0, +2, -2, +3, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*    2*/ -6.353000000000e+01,
      +2.000000000000e-02,
      {+0, +0, +0, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*    1*/ -2.640730000000e+03,
      +3.900000000000e-01,
      {+0, +0, +0, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0,
       +0}}}}; /* C0Series */

/* Expansion for s(t) series; j=1, Number of terms = 3
 * Extracted from Table 5.2c
 * Note that we are storing elements in reverse order (i.e. from i=3 to 1),
 * to traverse the smallest terms first.
 */
constexpr const std::array<CioSeriesData, 3> C1Series = {
    {{/*   36*/ +0.000000000000e+00,
      +4.800000000000e-01,
      {+0, +0, +2, -2, +3, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   35*/ +1.730000000000e+00,
      -3.000000000000e-02,
      {+0, +0, +0, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   34*/ -7.000000000000e-02,
      +3.570000000000e+00,
      {+0, +0, +0, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0,
       +0}}}}; /* C1Series */

/* Expansion for s(t) series; j=2, Number of terms = 25
 * Extracted from Table 5.2c
 * Note that we are storing elements in reverse order (i.e. from i=25 to 1),
 * to traverse the smallest terms first.
 */
constexpr const std::array<CioSeriesData, 25> C2Series = {
    {{/*   61*/ -1.100000000000e-01,
      +0.000000000000e+00,
      {+0, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   60*/ -1.200000000000e-01,
      +0.000000000000e+00,
      {+1, +0, +2, -2, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   59*/ -1.300000000000e-01,
      +0.000000000000e+00,
      {+2, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   58*/ +1.300000000000e-01,
      +0.000000000000e+00,
      {+2, +0, +2, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   57*/ +1.700000000000e-01,
      +0.000000000000e+00,
      {+0, +0, +2, +2, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   56*/ +2.000000000000e-01,
      +0.000000000000e+00,
      {+2, +0, -2, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   55*/ -2.100000000000e-01,
      +0.000000000000e+00,
      {+2, +0, +0, -2, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   54*/ +2.200000000000e-01,
      +0.000000000000e+00,
      {+1, +0, +2, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   53*/ -2.500000000000e-01,
      +0.000000000000e+00,
      {+1, +0, +0, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   52*/ -2.600000000000e-01,
      +0.000000000000e+00,
      {+1, +0, -2, -2, -2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   51*/ -2.700000000000e-01,
      +0.000000000000e+00,
      {+1, +0, +0, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   50*/ -2.700000000000e-01,
      +0.000000000000e+00,
      {+0, +0, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   49*/ +5.300000000000e-01,
      +0.000000000000e+00,
      {+1, +0, -2, +0, -2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   48*/ -5.500000000000e-01,
      +0.000000000000e+00,
      {+0, +0, +2, -2, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   47*/ +6.800000000000e-01,
      +0.000000000000e+00,
      {+1, +0, +0, -2, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   46*/ +9.300000000000e-01,
      +0.000000000000e+00,
      {+0, +1, -2, +2, -2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   45*/ +1.300000000000e+00,
      +0.000000000000e+00,
      {+1, +0, +2, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   44*/ +1.670000000000e+00,
      +0.000000000000e+00,
      {+0, +0, +2, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   43*/ +2.230000000000e+00,
      +0.000000000000e+00,
      {+0, +1, +2, -2, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   42*/ -3.070000000000e+00,
      +0.000000000000e+00,
      {+1, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   41*/ -6.380000000000e+00,
      -5.000000000000e-02,
      {+0, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   40*/ -8.850000000000e+00,
      +1.000000000000e-02,
      {+0, +0, +0, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   39*/ +9.840000000000e+00,
      -1.000000000000e-02,
      {+0, +0, +2, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   38*/ +5.691000000000e+01,
      +6.000000000000e-02,
      {+0, +0, +2, -2, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   37*/ +7.435200000000e+02,
      -1.700000000000e-01,
      {+0, +0, +0, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0,
       +0}}}}; /* C2Series */

/* Expansion for s(t) series; j=3, Number of terms = 4
 * Extracted from Table 5.2c
 * Note that we are storing elements in reverse order (i.e. from i=4 to 1),
 * to traverse the smallest terms first.
 */
constexpr const std::array<CioSeriesData, 4> C3Series = {
    {{/*   65*/ +0.000000000000e+00,
      +2.300000000000e-01,
      {+0, +0, +0, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   64*/ -1.000000000000e-02,
      -2.500000000000e-01,
      {+0, +0, +2, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   63*/ -3.000000000000e-02,
      -1.460000000000e+00,
      {+0, +0, +2, -2, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
     {/*   62*/ +3.000000000000e-01,
      -2.342000000000e+01,
      {+0, +0, +0, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0,
       +0}}}}; /* C3Series */

constexpr const std::array<CioSeriesData, 1> C4Series = {{
    {/*   66*/ -2.600000000000e-01,
     -1.000000000000e-02,
     {+0, +0, +0, +0, +1, +0, +0, +0, +0, +0, +0, +0, +0, +0}},
}}; /* C4Series */
} /* unnamed namespace */

double dso::detail::s06(const double *const fargs, double t, double x,
                        double y) noexcept {
  /* Compute the development for s + X*Y/2 */
  /* sum for j=4 */
  double c4(0e0);
  for (const auto &s : C4Series) {
    /* compute ARGUMENT and its trigs */
    const double arg = s.arg(fargs);
    const double ca = std::cos(arg);
    const double sa = std::sin(arg);
    /* add contribution from this frequency (i.e. i) */
    c4 += s.ccj * ca + s.csj * sa;
  }

  /* sum for j=3 */
  double c3(0e0);
  for (const auto &s : C3Series) {
    /* compute ARGUMENT and its trigs */
    const double arg = s.arg(fargs);
    const double ca = std::cos(arg);
    const double sa = std::sin(arg);
    /* add contribution from this frequency (i.e. i) */
    c3 += s.ccj * ca + s.csj * sa;
  }

  /* sum for j=2 */
  double c2(0e0);
  for (const auto &s : C2Series) {
    /* compute ARGUMENT and its trigs */
    const double arg = s.arg(fargs);
    const double ca = std::cos(arg);
    const double sa = std::sin(arg);
    /* add contribution from this frequency (i.e. i) */
    c2 += s.ccj * ca + s.csj * sa;
  }

  /* sum for j=1 */
  double c1(0e0);
  for (const auto &s : C1Series) {
    /* compute ARGUMENT and its trigs */
    const double arg = s.arg(fargs);
    const double ca = std::cos(arg);
    const double sa = std::sin(arg);
    /* add contribution from this frequency (i.e. i) */
    c1 += s.ccj * ca + s.csj * sa;
  }

  /* sum for j=0 */
  double c0(0e0);
  for (const auto &s : C0Series) {
    /* compute ARGUMENT and its trigs */
    const double arg = s.arg(fargs);
    const double ca = std::cos(arg);
    const double sa = std::sin(arg);
    /* add contribution from this frequency (i.e. i) */
    c0 += s.ccj * ca + s.csj * sa;
  }

  /* Add polynomial terms in [microarcseconds] */
  constexpr const double CPolyCoeffs[] = {94e0,        3808.65e0, -122.68e0,
                                          -72574.11e0, 27.98e0,   15.62e0};
  const double cpoly =
      CPolyCoeffs[0] +
      (CPolyCoeffs[1] +
       (CPolyCoeffs[2] +
        (CPolyCoeffs[3] + (CPolyCoeffs[4] + CPolyCoeffs[5] * t) * t) * t) *
           t) *
          t;

  /* accumulate in [microarcseconds] */
  const auto cseries = c0 + (c1 + (c2 + (c3 + c4 * t) * t) * t) * t;

  /* s + XY/2 angle in radians */
  return dso::sec2rad((cseries + cpoly) * 1e-6) - x * y / 2e0;
}

double dso::s06(const dso::MjdEpoch &tt, double x, double y) noexcept {
  /* Interval between fundamental date J2000.0 and given date. */
  const double t = tt.jcenturies_sinceJ2000();

  /* Fundamental arguments (IERS 2003) */
  const double fa[14] = {
      dso::iers2010::fal03(t),  dso::iers2010::falp03(t),
      dso::iers2010::faf03(t),  dso::iers2010::fad03(t),
      dso::iers2010::faom03(t), dso::iers2010::fame03(t),
      dso::iers2010::fave03(t), dso::iers2010::fae03(t),
      dso::iers2010::fama03(t), dso::iers2010::faju03(t),
      dso::iers2010::fasa03(t), dso::iers2010::faur03(t),
      dso::iers2010::fane03(t), dso::iers2010::fapa03(t),
  };
  return dso::detail::s06(fa, t, x, y);
}

double dso::s06(const dso::MjdEpoch &tt, double x, double y,
                const double *const fargs14) noexcept {
  return dso::detail::s06(fargs14, tt.jcenturies_sinceJ2000(), x, y);
}
