#include "eop.hpp"
#include "fundarg.hpp"

void dso::EopSeries::regularize_impl(int addremove) noexcept {
  double fargs[5];
  double ut1_cor, lod_cor, dummy;

  /* are we adding or removing the contribution ? */
  const double f = (addremove > 0) ? 1e-6 : -1e-6;

  for (auto &entry : mvec) {
    /* compute fundamental arguments */
    dso::fundarg(entry.t(), fargs);
    /* zonal tide corrections */
    dso::deop_zonal_tide(fargs, ut1_cor, lod_cor, dummy);
    /* apply (note: microseconds to seconds) */
    entry.dut() += (ut1_cor * f);
    entry.lod() += (lod_cor * f);
  }

  return;
}
