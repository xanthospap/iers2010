#include "iau.hpp"
#include "iersc.hpp"

// ----------------------------------------- //
// The series for the EE complementary terms //
// ----------------------------------------- //

typedef struct {
  int nfa[8];  // coefficients of l,l',F,D,Om,LVe,LE,pA
  double s, c; // sine and cosine coefficients
} TERM;

// Terms of order t^0
constexpr TERM e0[] = {

    /* 1-10 */
    {{0, 0, 0, 0, 1, 0, 0, 0}, 2640.96e-6, -0.39e-6},
    {{0, 0, 0, 0, 2, 0, 0, 0}, 63.52e-6, -0.02e-6},
    {{0, 0, 2, -2, 3, 0, 0, 0}, 11.75e-6, 0.01e-6},
    {{0, 0, 2, -2, 1, 0, 0, 0}, 11.21e-6, 0.01e-6},
    {{0, 0, 2, -2, 2, 0, 0, 0}, -4.55e-6, 0.00e-6},
    {{0, 0, 2, 0, 3, 0, 0, 0}, 2.02e-6, 0.00e-6},
    {{0, 0, 2, 0, 1, 0, 0, 0}, 1.98e-6, 0.00e-6},
    {{0, 0, 0, 0, 3, 0, 0, 0}, -1.72e-6, 0.00e-6},
    {{0, 1, 0, 0, 1, 0, 0, 0}, -1.41e-6, -0.01e-6},
    {{0, 1, 0, 0, -1, 0, 0, 0}, -1.26e-6, -0.01e-6},

    /* 11-20 */
    {{1, 0, 0, 0, -1, 0, 0, 0}, -0.63e-6, 0.00e-6},
    {{1, 0, 0, 0, 1, 0, 0, 0}, -0.63e-6, 0.00e-6},
    {{0, 1, 2, -2, 3, 0, 0, 0}, 0.46e-6, 0.00e-6},
    {{0, 1, 2, -2, 1, 0, 0, 0}, 0.45e-6, 0.00e-6},
    {{0, 0, 4, -4, 4, 0, 0, 0}, 0.36e-6, 0.00e-6},
    {{0, 0, 1, -1, 1, -8, 12, 0}, -0.24e-6, -0.12e-6},
    {{0, 0, 2, 0, 0, 0, 0, 0}, 0.32e-6, 0.00e-6},
    {{0, 0, 2, 0, 2, 0, 0, 0}, 0.28e-6, 0.00e-6},
    {{1, 0, 2, 0, 3, 0, 0, 0}, 0.27e-6, 0.00e-6},
    {{1, 0, 2, 0, 1, 0, 0, 0}, 0.26e-6, 0.00e-6},

    /* 21-30 */
    {{0, 0, 2, -2, 0, 0, 0, 0}, -0.21e-6, 0.00e-6},
    {{0, 1, -2, 2, -3, 0, 0, 0}, 0.19e-6, 0.00e-6},
    {{0, 1, -2, 2, -1, 0, 0, 0}, 0.18e-6, 0.00e-6},
    {{0, 0, 0, 0, 0, 8, -13, -1}, -0.10e-6, 0.05e-6},
    {{0, 0, 0, 2, 0, 0, 0, 0}, 0.15e-6, 0.00e-6},
    {{2, 0, -2, 0, -1, 0, 0, 0}, -0.14e-6, 0.00e-6},
    {{1, 0, 0, -2, 1, 0, 0, 0}, 0.14e-6, 0.00e-6},
    {{0, 1, 2, -2, 2, 0, 0, 0}, -0.14e-6, 0.00e-6},
    {{1, 0, 0, -2, -1, 0, 0, 0}, 0.14e-6, 0.00e-6},
    {{0, 0, 4, -2, 4, 0, 0, 0}, 0.13e-6, 0.00e-6},

    /* 31-33 */
    {{0, 0, 2, -2, 4, 0, 0, 0}, -0.11e-6, 0.00e-6},
    {{1, 0, -2, 0, -3, 0, 0, 0}, 0.11e-6, 0.00e-6},
    {{1, 0, -2, 0, -1, 0, 0, 0}, 0.11e-6, 0.00e-6}};

// Terms of order t^1
constexpr TERM e1[] = {{{0, 0, 0, 0, 1, 0, 0, 0}, -0.87e-6, 0.00e-6}};

// Number of terms in the series
constexpr int NE0 = (int)(sizeof e0 / sizeof(TERM));
constexpr int NE1 = (int)(sizeof e1 / sizeof(TERM));

double iers2010::sofa::eect00(double date1, double date2) noexcept {

  // Interval between fundamental epoch J2000.0 and current date (JC).
  // Time since J2000.0, in Julian centuries
  const double t = ((date1 - iers2010::DJ00) + date2) / iers2010::DJC;

  // Fundamental Arguments (from IERS Conventions 2003)
  double fa[] = {// Mean anomaly of the Moon.
                 fal03(t),
                 // Mean anomaly of the Sun.
                 falp03(t),
                 // Mean longitude of the Moon minus that of the ascending node.
                 faf03(t),
                 // Mean elongation of the Moon from the Sun.
                 fad03(t),
                 // Mean longitude of the ascending node of the Moon.
                 faom03(t),
                 // Mean longitude of Venus.
                 fave03(t),
                 // Mean longitude of Earth.
                 fae03(t),
                 // General precession in longitude.
                 fapa03(t)};

  // Evaluate the EE complementary terms.
  double s0 = 0e0;
  double s1 = 0e0;

  for (int i = NE0 - 1; i >= 0; i--) {
    double a = 0e0;
    for (int j = 0; j < 8; j++) {
      a += (double)(e0[i].nfa[j]) * fa[j];
    }
    s0 += e0[i].s * std::sin(a) + e0[i].c * std::cos(a);
  }

  for (int i = NE1 - 1; i >= 0; i--) {
    double a = 0e0;
    for (int j = 0; j < 8; j++) {
      a += (double)(e1[i].nfa[j]) * fa[j];
    }
    s1 += e1[i].s * std::sin(a) + e1[i].c * std::cos(a);
  }

  return (s0 + s1 * t) * iers2010::DAS2R;
}
