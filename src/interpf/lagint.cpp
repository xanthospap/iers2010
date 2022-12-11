#include "iers2010.hpp"

int iers2010::interp::lagint(const double *x, const double *y, int n,
                             double xint, double &yout, int &idx,
                             int order) noexcept {
  yout = 0e0;
  
  const int window = order + 1;
  if ( n<window || window%2 ) return 1;

  // k is the index for which xint >= x[k] and xint < x[k+1]
  int k = idx;

  // if k not provided (aka k<0), find the interval and k value
  if (k < 0 || !(xint >= x[k] && xint < x[k + 1])) {
    for (int i = 0; i < n - 1; i++) {
      if (xint >= x[i] && xint < x[i + 1]) {
        k = i;
        break;
      }
    }
  }

  // we are going to use window points for the interpolation
  // at least window/2 points before and window/2 after should be available
  const int npts = window / 2;
  int status = 0;
  status += !(k>=npts-1) + !(k<n-npts);
  if (status) return status;


  for (int m = k + 1 - npts; m <= k + npts; m++) {
    double term = y[m];
    for (int j = k + 1 - npts; j < m; j++)
      term *= (xint - x[j]) / (x[m] - x[j]);
    for (int j = m + 1; j <= k + npts; j++)
      term *= (xint - x[j]) / (x[m] - x[j]);
    yout += term;
  }

  idx = k;
  return 0;
}
