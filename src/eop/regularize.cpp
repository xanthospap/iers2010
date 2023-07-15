#include "eop.hpp"
#include <algorithm>

void dso::EopLookUpTable::regularizedCorrections(std::vector<double> &dUt1R,
                                                 std::vector<double> &dLodR,
                                                 bool tInUtc) const noexcept {
  /* clear and reserve */
  dUt1R.clear();
  dUt1R.reserve(tvec.size());
  dLodR.clear();
  dLodR.reserve(tvec.size());

  /* hold fundamental arguments plus GMST */
  double fargs[6];

  /* iterator to the ΔUT1 values (matching dates) */
  auto dut1it = dut1.begin();

  /* compute zonal tides for ΔUT1 and ΔLOD */
  double ut1r, lodr, om;
  for (auto epoch = tvec.cbegin(); epoch != tvec.cend(); ++epoch) {
    /* UTC to TT for current (tabulated date) */
    const auto tt = (tInUtc) ? (epoch->utc2tt()) : (*epoch);
    /* fundamental arguments */
    iers2010::fundarg(tt, fargs);
    /* Compute effect of zonal harmonics using RG_ZONT2 */
    iers2010::rg_zont2(fargs, ut1r, lodr, om);
    dUt1R.push_back(ut1r);
    dLodR.push_back(lodr);
    /* augment ΔUT1 iterator */
    ++dut1it;
  }

  /* all done */
  return;
}

void dso::EopLookUpTable::regularize(bool tInUtc) noexcept {
  std::vector<double> dUt1R, dLodR;
  dUt1R.reserve(tvec.size());
  dLodR.reserve(tvec.size());

  /* compute zonal tide corrections */
  this->regularizedCorrections(dUt1R, dLodR, tInUtc);

  /* apply corrections to get "regularized" values */
  std::transform(
      dut1.cbegin(), dut1.cend(), dUt1R.cbegin(), dut1.begin(),
      [](const double dref, const double dcor) { return dref - dcor; });
  std::transform(
      lod.cbegin(), lod.cend(), dLodR.cbegin(), lod.begin(),
      [](const double dref, const double dcor) { return dref - dcor; });

  return;
}
