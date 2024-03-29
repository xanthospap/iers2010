#include "iau.hpp"

namespace {

/* The series for s+XY/2 */
typedef struct {
  int nfa[8];  /* coefficients of l,l',F,D,Om,LVe,LE,pA */
  double s, c; /* sin and cos coefficients */
} TERM;

/* Polynomial coefficients */
constexpr const double sp[] = {
    /* 1-6 */
    94.00e-6, 3808.35e-6, -119.94e-6, -72574.09e-6, 27.70e-6, 15.61e-6};

/* Terms of order t^0 */
constexpr const TERM s0[] = {

    /* 1-10 */
    {{0, 0, 0, 0, 1, 0, 0, 0}, -2640.73e-6, 0.39e-6},
    {{0, 0, 0, 0, 2, 0, 0, 0}, -63.53e-6, 0.02e-6},
    {{0, 0, 2, -2, 3, 0, 0, 0}, -11.75e-6, -0.01e-6},
    {{0, 0, 2, -2, 1, 0, 0, 0}, -11.21e-6, -0.01e-6},
    {{0, 0, 2, -2, 2, 0, 0, 0}, 4.57e-6, 0.00e-6},
    {{0, 0, 2, 0, 3, 0, 0, 0}, -2.02e-6, 0.00e-6},
    {{0, 0, 2, 0, 1, 0, 0, 0}, -1.98e-6, 0.00e-6},
    {{0, 0, 0, 0, 3, 0, 0, 0}, 1.72e-6, 0.00e-6},
    {{0, 1, 0, 0, 1, 0, 0, 0}, 1.41e-6, 0.01e-6},
    {{0, 1, 0, 0, -1, 0, 0, 0}, 1.26e-6, 0.01e-6},

    /* 11-20 */
    {{1, 0, 0, 0, -1, 0, 0, 0}, 0.63e-6, 0.00e-6},
    {{1, 0, 0, 0, 1, 0, 0, 0}, 0.63e-6, 0.00e-6},
    {{0, 1, 2, -2, 3, 0, 0, 0}, -0.46e-6, 0.00e-6},
    {{0, 1, 2, -2, 1, 0, 0, 0}, -0.45e-6, 0.00e-6},
    {{0, 0, 4, -4, 4, 0, 0, 0}, -0.36e-6, 0.00e-6},
    {{0, 0, 1, -1, 1, -8, 12, 0}, 0.24e-6, 0.12e-6},
    {{0, 0, 2, 0, 0, 0, 0, 0}, -0.32e-6, 0.00e-6},
    {{0, 0, 2, 0, 2, 0, 0, 0}, -0.28e-6, 0.00e-6},
    {{1, 0, 2, 0, 3, 0, 0, 0}, -0.27e-6, 0.00e-6},
    {{1, 0, 2, 0, 1, 0, 0, 0}, -0.26e-6, 0.00e-6},

    /* 21-30 */
    {{0, 0, 2, -2, 0, 0, 0, 0}, 0.21e-6, 0.00e-6},
    {{0, 1, -2, 2, -3, 0, 0, 0}, -0.19e-6, 0.00e-6},
    {{0, 1, -2, 2, -1, 0, 0, 0}, -0.18e-6, 0.00e-6},
    {{0, 0, 0, 0, 0, 8, -13, -1}, 0.10e-6, -0.05e-6},
    {{0, 0, 0, 2, 0, 0, 0, 0}, -0.15e-6, 0.00e-6},
    {{2, 0, -2, 0, -1, 0, 0, 0}, 0.14e-6, 0.00e-6},
    {{0, 1, 2, -2, 2, 0, 0, 0}, 0.14e-6, 0.00e-6},
    {{1, 0, 0, -2, 1, 0, 0, 0}, -0.14e-6, 0.00e-6},
    {{1, 0, 0, -2, -1, 0, 0, 0}, -0.14e-6, 0.00e-6},
    {{0, 0, 4, -2, 4, 0, 0, 0}, -0.13e-6, 0.00e-6},

    /* 31-33 */
    {{0, 0, 2, -2, 4, 0, 0, 0}, 0.11e-6, 0.00e-6},
    {{1, 0, -2, 0, -3, 0, 0, 0}, -0.11e-6, 0.00e-6},
    {{1, 0, -2, 0, -1, 0, 0, 0}, -0.11e-6, 0.00e-6}};

/* Terms of order t^1 */
constexpr const TERM s1[] = {

    /* 1-3 */
    {{0, 0, 0, 0, 2, 0, 0, 0}, -0.07e-6, 3.57e-6},
    {{0, 0, 0, 0, 1, 0, 0, 0}, 1.71e-6, -0.03e-6},
    {{0, 0, 2, -2, 3, 0, 0, 0}, 0.00e-6, 0.48e-6}};

/* Terms of order t^2 */
constexpr const TERM s2[] = {

    /* 1-10 */
    {{0, 0, 0, 0, 1, 0, 0, 0}, 743.53e-6, -0.17e-6},
    {{0, 0, 2, -2, 2, 0, 0, 0}, 56.91e-6, 0.06e-6},
    {{0, 0, 2, 0, 2, 0, 0, 0}, 9.84e-6, -0.01e-6},
    {{0, 0, 0, 0, 2, 0, 0, 0}, -8.85e-6, 0.01e-6},
    {{0, 1, 0, 0, 0, 0, 0, 0}, -6.38e-6, -0.05e-6},
    {{1, 0, 0, 0, 0, 0, 0, 0}, -3.07e-6, 0.00e-6},
    {{0, 1, 2, -2, 2, 0, 0, 0}, 2.23e-6, 0.00e-6},
    {{0, 0, 2, 0, 1, 0, 0, 0}, 1.67e-6, 0.00e-6},
    {{1, 0, 2, 0, 2, 0, 0, 0}, 1.30e-6, 0.00e-6},
    {{0, 1, -2, 2, -2, 0, 0, 0}, 0.93e-6, 0.00e-6},

    /* 11-20 */
    {{1, 0, 0, -2, 0, 0, 0, 0}, 0.68e-6, 0.00e-6},
    {{0, 0, 2, -2, 1, 0, 0, 0}, -0.55e-6, 0.00e-6},
    {{1, 0, -2, 0, -2, 0, 0, 0}, 0.53e-6, 0.00e-6},
    {{0, 0, 0, 2, 0, 0, 0, 0}, -0.27e-6, 0.00e-6},
    {{1, 0, 0, 0, 1, 0, 0, 0}, -0.27e-6, 0.00e-6},
    {{1, 0, -2, -2, -2, 0, 0, 0}, -0.26e-6, 0.00e-6},
    {{1, 0, 0, 0, -1, 0, 0, 0}, -0.25e-6, 0.00e-6},
    {{1, 0, 2, 0, 1, 0, 0, 0}, 0.22e-6, 0.00e-6},
    {{2, 0, 0, -2, 0, 0, 0, 0}, -0.21e-6, 0.00e-6},
    {{2, 0, -2, 0, -1, 0, 0, 0}, 0.20e-6, 0.00e-6},

    /* 21-25 */
    {{0, 0, 2, 2, 2, 0, 0, 0}, 0.17e-6, 0.00e-6},
    {{2, 0, 2, 0, 2, 0, 0, 0}, 0.13e-6, 0.00e-6},
    {{2, 0, 0, 0, 0, 0, 0, 0}, -0.13e-6, 0.00e-6},
    {{1, 0, 2, -2, 2, 0, 0, 0}, -0.12e-6, 0.00e-6},
    {{0, 0, 2, 0, 0, 0, 0, 0}, -0.11e-6, 0.00e-6}};

/* Terms of order t^3 */
constexpr const TERM s3[] = {

    /* 1-4 */
    {{0, 0, 0, 0, 1, 0, 0, 0}, 0.30e-6, -23.51e-6},
    {{0, 0, 2, -2, 2, 0, 0, 0}, -0.03e-6, -1.39e-6},
    {{0, 0, 2, 0, 2, 0, 0, 0}, -0.01e-6, -0.24e-6},
    {{0, 0, 0, 0, 2, 0, 0, 0}, 0.00e-6, 0.22e-6}};

/* Terms of order t^4 */
constexpr const TERM s4[] = {

    /* 1-1 */
    {{0, 0, 0, 0, 1, 0, 0, 0}, -0.26e-6, -0.01e-6}};

/* Number of terms in the series */
constexpr int NS0 = (int)(sizeof s0 / sizeof(TERM));
constexpr int NS1 = (int)(sizeof s1 / sizeof(TERM));
constexpr int NS2 = (int)(sizeof s2 / sizeof(TERM));
constexpr int NS3 = (int)(sizeof s3 / sizeof(TERM));
constexpr int NS4 = (int)(sizeof s4 / sizeof(TERM));
} /* unnamed namespace */

double iers2010::sofa::s00(const dso::TwoPartDate &mjd_tt, double x,
                           double y) noexcept {
  /* Interval between fundamental epoch J2000.0 and current date (JC). */
  const double t = mjd_tt.jcenturies_sinceJ2000();

  /* Fundamental Arguments (from IERS Conventions 2003) */
  const double fa[8] = {
      /* Mean anomaly of the Moon. */
      iers2010::fal03(t),
      /* Mean anomaly of the Sun. */
      iers2010::falp03(t),
      /* Mean longitude of the Moon minus that of the ascending node. */
      iers2010::faf03(t),
      /* Mean elongation of the Moon from the Sun. */
      iers2010::fad03(t),
      /* Mean longitude of the ascending node of the Moon. */
      iers2010::faom03(t),
      /* Mean longitude of Venus. */
      iers2010::fave03(t),
      /* Mean longitude of Earth. */
      iers2010::fae03(t),
      /* General precession in longitude. */
      iers2010::fapa03(t)};

  /* Evaluate s. */
  double w0 = sp[0];
  double w1 = sp[1];
  double w2 = sp[2];
  double w3 = sp[3];
  double w4 = sp[4];
  double w5 = sp[5];

  for (int i = NS0 - 1; i >= 0; i--) {
    double a = 0e0;
    for (int j = 0; j < 8; j++) {
      a += s0[i].nfa[j] * fa[j];
    }
    const double sa = std::sin(a);
    const double ca = std::cos(a);
    w0 += s0[i].s * sa + s0[i].c * ca;
  }

  for (int i = NS1 - 1; i >= 0; i--) {
    double a = 0e0;
    for (int j = 0; j < 8; j++) {
      a += s1[i].nfa[j] * fa[j];
    }
    const double sa = std::sin(a);
    const double ca = std::cos(a);
    w1 += s1[i].s * sa + s1[i].c * ca;
  }

  for (int i = NS2 - 1; i >= 0; i--) {
    double a = 0e0;
    for (int j = 0; j < 8; j++) {
      a += s2[i].nfa[j] * fa[j];
    }
    const double sa = std::sin(a);
    const double ca = std::cos(a);
    w2 += s2[i].s * sa + s2[i].c * ca;
  }

  for (int i = NS3 - 1; i >= 0; i--) {
    double a = 0e0;
    for (int j = 0; j < 8; j++) {
      a += s3[i].nfa[j] * fa[j];
    }
    const double sa = std::sin(a);
    const double ca = std::cos(a);
    w3 += s3[i].s * sa + s3[i].c * ca;
  }

  for (int i = NS4 - 1; i >= 0; i--) {
    double a = 0e0;
    for (int j = 0; j < 8; j++) {
      a += s4[i].nfa[j] * fa[j];
    }
    const double sa = std::sin(a);
    const double ca = std::cos(a);
    w4 += s4[i].s * sa + s4[i].c * ca;
  }

  const double s =
      dso::sec2rad(w0 + (w1 + (w2 + (w3 + (w4 + w5 * t) * t) * t) * t) * t) -
      x * y / 2e0;

  return s;
}
