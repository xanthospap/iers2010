#include "iers2010.hpp"
#include <cstdio>

int iers2010::interp_pole(const double *mjd, const double *x, const double *y,
                          const double *ut1, int n, double rjd, double &xint,
                          double &yint, double &ut1int) noexcept {
  // quick return in case the passed in date is out of bounds ...
  if (rjd<*mjd || rjd>=mjd[n-1]) return 1;

  int status = 0;
  int index = -1;
  status = iers2010::interp::lagint(mjd, x, n, rjd, xint, index);
  status = iers2010::interp::lagint(mjd, y, n, rjd, yint, index);
  status = iers2010::interp::lagint(mjd, ut1, n, rjd, ut1int, index);

  // corrections
  double corx, cory, corut1, corlod;

  // compute julian centuries
  const double tjc = (rjd - 51544.5e0)/36525.0e0;

  // oceanic effect
  iers2010::interp::pmut1_oceans(tjc, corx, cory, corut1, corlod);
  xint += corx;
  yint += cory;
  ut1int += corut1;

  // lunisolar effect
  iers2010::interp::pm_gravi(tjc, corx, cory);
  xint += corx;
  yint += cory;

  return status;
}
