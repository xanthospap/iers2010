#include "iers2010.hpp"

int iers2010::interp::lagint(const double *x, const double *y, int n,
                             double xint, double &yout, int &idx) noexcept {
  int status = 0;
  yout = 0e0;
  int k = idx;
  if (k < 0) {
    for (int i = 0; i < n - 1; i++) {
      if (xint >= x[i] && xint < x[i + 1]) {
        k = i;
        break;
      }
    }
  }

  if (k < 1) {
    k = 1;
    status = -1;
  }

  if (k > n - 3) {
    k = n - 3;
    status = -2;
  }

  for (int m = k - 1; m <= k + 2; m++) {
    double term = y[m];
    for (int j = k - 1; j <= k + 2; j++) {
      if (m != j){
        term *= (xint - x[j]) / (x[m] - x[j]);
      }
    }
    yout += term;
  }

  idx = k;
  return status;
}
