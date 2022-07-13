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
  // printf("\tNote: Interpolated values: %+10.2f %+10.2f %+10.2f\n", xint, yint, ut1int);

  // corrections
  double corx, cory, corut1, corlod;

  // compute julian centuries
  const double tjc = (rjd - 51544.5e0)/36525.0e0;

  // oceanic effect (results in [μas])
  iers2010::interp::pmut1_oceans(tjc, corx, cory, corut1, corlod);
  xint += (corx * 1e-3);
  yint += (cory * 1e-3);
  ut1int += (corut1 * 1e-3);
  // printf("\tNote: Oceanic Effect     : %+10.3f %+10.3f %+10.3f\n", corx*1e-3, cory*1e-3, corut1*1e-3);

  // lunisolar effect (results in [μas])
  iers2010::interp::pm_gravi(tjc, corx, cory);
  xint += (corx * 1e-3);
  yint += (cory * 1e-3);
  // printf("\tNote: LuniSolar Effect   : %+10.3f %+10.3f %+10.3f\n", corx*1e-3, cory*1e-3, 0e0);

  return status;
}
