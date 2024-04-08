#include "eop.hpp"

void dso::EopSeries::regularize() noexcept {
  double fargs[5];
  double ut1_cor, lod_cor, dummy;
  for (auto &entry : mvec) {
    /* compute fundamental arguments */
    dso::fundarg(entry.t(), fargs);
    /* zonal tide corrections */
    dso::deop_zonal_tide(fargs, ut1_cor, lod_cor, dummy);
    /* apply (note: microseconds to seconds) */
    entry.dut() += (ut1_cor * 1e-6);
    entry.lod() += (lod_cor * 1e-6);
  }

  return;
}
