#include "gravity.hpp"
#include <cassert>
#include <cmath>

namespace {
inline double anm(int n, int m) noexcept {
  return std::sqrt( (double)((2*n-1)*(2*n+1)) / (double)((n-m)*(n+m)) );
}

inline double bnm(int n, int m) noexcept {
  long b1 = (2*n+1)*(n+m-1)*(n-m-1);
  long b2 = (n-m)*(n+m)*(2*n-3);
  return std::sqrt( (double)b1 / (double)b2 );
}
} /* unnamed namespace */

int dso::normalised_associated_legendre_functions(
    double theta, int max_degree, int max_order,
    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise>
        &P) noexcept {
  /* check degree and order */
  assert(max_order <= max_degree);

  /* trigs for expansion */
  const double u = std::cos(theta);
  const double t = std::sin(theta);

  /* first two terms */
  P(0, 0) = 1e0;
  P(1, 1) = std::sqrt(3e0) * u;

  /* loop through columns, m=[0,max_order) */
  for (int m = 0; m < max_order; m++) {
    /* first value from sectorial diagonal, i.e. P(m+1,m) has no n-2 term */
    int n = m + 1;
    double *__restrict__ pmm = P.column(m);
    // P(n,m) = anm(n,m) * t * P(m,m);
    pmm[1] = anm(n, m) * t * pmm[0];
    /* all other terms down for col=m, i.e. n=[m+2,max_degree] */
    for (n = m + 2; n <= max_degree; n++) {
      //printf("\tP(%d,%d) = %+15.9f*%.6f %+15.9f*%.6f\n", n,m,anm(n,m), bnm(n,m),P(n-1,m),P(n-2,m));
      //P(n,m) = anm(n,m) * t * P(n-1,m) - bnm(n,m) * P(n-2,m);
      int k = n - m;
      pmm[k] = anm(n, m) * t * pmm[k - 1] - bnm(n, m) * pmm[k - 2];
    }
    /* next sectorial term P(m+1,m+1) */
    if (m>0) 
    // P(m+1,m+1) = u * std::sqrt( (2e0*m+3)/(2e0*m+2) ) * P(m,m);
    P.column(m + 1)[0] = u * std::sqrt((2e0 * m + 3) / (2e0 * m+2)) * pmm[0];
  }

  /* last column, m=max_order */
  int m = max_order;
  double *__restrict__ pmm = P.column(m);
  pmm[1] = anm(m + 1, m) * t * pmm[0];
  for (int n = m + 2; n <= max_degree; n++) {
    int k = n - m;
    pmm[k] = anm(n, m) * t * pmm[k - 1] - bnm(n, m) * pmm[k - 2];
  }

  return 0;
}

//int dso::modified_normalised_associated_legendre_functions(
//    double theta, int max_degree, int max_order,
//    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> &P,
//    dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise>
//        &Pdot) noexcept {
//  /* check degree and order */
//  assert(max_order <= max_degree);
//
//  /* trigs for expansion */
//  const double u = std::cos(theta);
//  const double t = std::sin(theta);
//
//  /* first two terms */
//  P(0, 0) = 0e0;
//  P(1, 1) = std::sqrt(3e0);
//  double um = 1e0; /* i.e. u^m */
//
//  /* loop through columns, m=[0,max_order) */
//  for (int m = 0; m < max_order; m++) {
//    /* first value from sectorial diagonal, i.e. P(m+1,m) has no n-2 term */
//    int n = m + 1;
//    double *__restrict__ pmm = P.column(m);
//    pmm[1] = anm(n, m) * t * pmm[0];
//    /* all other terms down for col=m, i.e. n=[m+2,max_degree] */
//    for (n = m + 2; n <= max_degree; n++) {
//      int k = n - m;
//      pmm[k] = anm(n, m) * t * pmm[k - 1] - bnm(n, m) * pmm[k - 2];
//    }
//    /* next sectorial term P(m+1,m+1) */
//    P.column(m + 1)[0] = std::sqrt((2e0 * m + 1) / (2e0 * m)) * pmm[0];
//  }
//
//  /* last column, m=max_order */
//  int m = max_order;
//  double *__restrict__ pmm = P.column(m);
//  pmm[1] = anm(m + 1, m) * t * pmm[0];
//  for (int n = m + 2; n <= max_degree; n++) {
//    int k = n - m;
//    pmm[k] = anm(n, m) * t * pmm[k - 1] - bnm(n, m) * pmm[k - 2];
//  }
//
//  return 0;
//}
