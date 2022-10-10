#include "iers2010.hpp"
#include <cstdio>

int iers2010::interp_pole(const double *mjd, const double *x, const double *y,
                          const double *ut1, int n, double rjd, double &xint,
                          double &yint, double &ut1int,
                          double &corlod) noexcept {
  // quick return in case the passed in date is out of bounds ...
  if (rjd<*mjd || rjd>=mjd[n-1]) {
    fprintf(stderr,
            "[ERROR] Requested interpolation at %.5f but EOP values span [%.5f "
            "to %.5f] (traceback: %s)\n",
            rjd, mjd[0], mjd[n - 1], __func__);
    return 1;
  }

  int status = 0;
  int index = -1;
  status += iers2010::interp::lagint(mjd, x, n, rjd, xint, index);
  status += iers2010::interp::lagint(mjd, y, n, rjd, yint, index);
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

  // compute julian centuries
  const double tjc = (rjd - 51544.5e0)/36525.0e0;

  // oceanic effect (results in [μas])
  iers2010::interp::pmut1_oceans(tjc, corx, cory, corut1, corlod);
  xint += (corx * 1e-3);
  yint += (cory * 1e-3);
  ut1int += (corut1 * 1e-3);
  corlod *= 1e-3;

  // lunisolar effect (results in [μas])
  iers2010::interp::pm_gravi(tjc, corx, cory);
  xint += (corx * 1e-3);
  yint += (cory * 1e-3);

  return status;
}
