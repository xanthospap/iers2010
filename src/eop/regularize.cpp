#include "eop.hpp"
#include <algorithm>

void dso::EopLookUpTable::regularizedCorrections(
    std::vector<double> &dUt1R, std::vector<double> &dLodR, bool tInUtc) const noexcept {
  // clear and reserve
  dUt1R.clear();
  dUt1R.reserve(t.size());
  dLodR.clear();
  dLodR.reserve(t.size());

  // compute zonal tides for ΔUT1 and ΔLOD
  double ut1r, lodr, om;
  for (auto epoch = t.cbegin(); epoch != t.cend(); ++epoch) {
    // UTC to TT for current (tabulated date)
    const auto tt = (tInUtc) ? (epoch->utc2tt()) : (*epoch);
    // Compute effect of zonal harmonics using RG_ZONT2
    iers2010::rg_zont2(tt, ut1r, lodr, om);
    dUt1R.push_back(ut1r);
    dLodR.push_back(lodr);
  }
  return;
}

void dso::EopLookUpTable::regularize(bool tInUtc) noexcept {
  std::vector<double> dUt1R, dLodR;
  dUt1R.reserve(t.size());
  dLodR.reserve(t.size());

  // compute zonal tide corrections
  this->regularizedCorrections(dUt1R, dLodR, tInUtc);

  // apply corrections to get "regularized" values
  std::transform(
      dut1.cbegin(), dut1.cend(), dUt1R.cbegin(), dut1.begin(),
      [](const double dref, const double dcor) { return dref - dcor; });
  std::transform(
      lod.cbegin(), lod.cend(), dLodR.cbegin(), lod.begin(),
      [](const double dref, const double dcor) { return dref - dcor; });

  return;
}
