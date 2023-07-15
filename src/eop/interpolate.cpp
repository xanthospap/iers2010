#include "eop.hpp"
#include <algorithm>

int dso::EopLookUpTable::interpolate(const dso::TwoPartDate &mjd,
                                     dso::EopRecord &eopr, int order,
                                     bool regularized) const noexcept {
  /* perform simple interpolation (lagrangian) */
  if (lagrange_interpolation(mjd, eopr, order))
    return 1;

  /* if these are "regularized" values for ΔUT1 and ΔLOD, add back the
   * zonal effects (for the requested epoch)
   */
  if (regularized) {
    double ut1r, lodr, om;
    iers2010::rg_zont2(mjd, ut1r, lodr, om);
    eopr.dut += ut1r;
    eopr.lod += lodr;
  }

  /* compute the diurnal and semi-diurnal variations in Earth Orientation
   * Parameters (x, y, UT1) from ocean tides.
   */
  double dxoc, dyoc, dut1oc; /* [μas] and [μsec] */
  iers2010::ortho_eop(mjd, dxoc, dyoc, dut1oc);

  /* compute fundamental arguments (and gmst+π) needed for pmsdnut2 and
   * utlibr
   */
  double fargs[6];
  iers2010::utils::fargs(mjd, eopr.dut, fargs);

  double dxlib, dylib; /* [μas] */
  iers2010::pmsdnut2(mjd, fargs, dxlib, dylib);

  double dut1lib, dlodlib;
  iers2010::utlibr(mjd, fargs, dut1lib, dlodlib);

  /* add corrections */
  eopr.xp += (dxoc + dxlib) * 1e-6;
  eopr.yp += (dyoc + dylib) * 1e-6;
  eopr.dut += (dut1oc + dut1lib) * 1e-6;
  eopr.lod += (dlodlib)*1e-6;

  return 0;
}
