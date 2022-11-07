#include "iers2010.hpp"
#include <cstdio>

int iers2010::interp_pole(const double *mjd, const double *x, const double *y,
                          const double *ut1, int n, double rjd, double &xint,
                          double &yint, double &ut1int,
                          double &corlod) noexcept {
  // perform lagrange interpolation (find index in first call and re-use)
  int status = 0;
  int index = -1;
  // x-pole ... (find index)
  status += iers2010::interp::lagint(mjd, x, n, rjd, xint, index);
  // y-pole ... (use index)
  status += iers2010::interp::lagint(mjd, y, n, rjd, yint, index);
  // dut ... (use index)
  status += iers2010::interp::lagint(mjd, ut1, n, rjd, ut1int, index);
  if (status) {
    fprintf(stderr,
            "[ERROR] Failed EOP interpolation; requested mjd %.5f (traceback: "
            "%s)\n",
            rjd, __func__);
    return 9;
  }

  // corrections
  double corx, cory, corut1;

  // compute julian centuries TT
  const double tjc = (rjd - 51544.5e0)/36525.0e0;

  // oceanic effect (results in [as] and [seconds])
  iers2010::interp::pmut1_oceans(tjc, corx, cory, corut1, corlod);
  xint += corx;
  yint += cory;
  ut1int += corut1;
  // corlod *= 1e-3;

  // lunisolar effect (results in [as])
  iers2010::interp::pm_gravi(tjc, corx, cory);
  xint += corx;
  yint += cory;

  return status;
}
