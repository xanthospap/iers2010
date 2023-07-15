#include "iau.hpp"
#include "iers2010.hpp"

int iers2010::utils::fargs(const dso::TwoPartDate &mjd_tt, double dut1,
                           double *fargs) noexcept {
  /* Evaluate the vector of the fundamental arguments
   * farg = [ GMST+pi, el, elp, f, d, om ] at t = fmjd
   */

  /* 1. Compute GMST using the IAU2006A model */
  const double gmst = iers2010::sofa::gmst06(mjd_tt.tt2ut1(dut1), mjd_tt);
  fargs[0] = gmst + iers2010::DPI;

  /* 2. Fundamental arguments -> GMST+pi, l, lp, f, d, om */
  fargs[1] = iers2010::fal03(mjd_tt);
  fargs[2] = iers2010::falp03(mjd_tt);
  fargs[3] = iers2010::faf03(mjd_tt);
  fargs[4] = iers2010::fad03(mjd_tt);
  fargs[5] = iers2010::faom03(mjd_tt);

  return 0;
}
