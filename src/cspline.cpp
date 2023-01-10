#include "cspline.hpp"

int dso::cspline_interp(double atx, int &index, const double *const x,
                   const double *const y, int n, const double *const y2,
                   double &yintrp) noexcept {
  if (index < 0) {
    // shit, we don;t know the right index for atx. Search for it!
    index = std::upper_bound(x, x + n, atx) - x;
    if (index >= n - 1) {
      return 99;
    }
  }
  assert(index < n - 1);
  const double h = x[index + 1] - x[index];
  const double a = (x[index + 1] - atx) / h;
  const double b = (atx - x[index]) / h;
  yintrp = a * y[index] + b * y[index + 1] +
           ((a * a * a - a) * y2[index] + (b * b * b - b) * y[index + 1]) *
               (h * h) / 6e0;
  return 0;
}

void dso::cspline_deriv(const double *const x, const double *const y, int n,
                   double *__restrict__ y2, double yp1, double ypn,
                   double *__restrict__ u) noexcept {

  // The lower boundary condition is set either to be "natural"
  // or else to have a speciﬁed ﬁrst derivative
  if (yp1 > .99e99) {
    y2[0] = u[0] = 0e0;
  } else {
    y2[0] = -.5e0;
    u[0] = (3e0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }

  for (int i = 1; i < n - 1; i++) {
    const double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    const double p = sig * y2[i - 1] + 2e0;
    y2[i] = (sig - 1e0) / p;
    u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
           (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    u[i] = (6e0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
  }

  double qn, un;
  if (ypn > .99e99) {
    qn = un = 0e0;
  } else {
    qn = .5e0;
    un = (3.0 / (x[n - 1] - x[n - 2])) *
         (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  }
  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1e0);

  // backsubstitution
  for (int k = n - 2; k >= 0; k--)
    y2[k] = y2[k] * y2[k + 1] + u[k];

  return;
}
