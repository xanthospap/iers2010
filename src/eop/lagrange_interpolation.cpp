#include "eop.hpp"
#include <algorithm>

namespace {
int lag_int(const dso::TwoPartDate &mjd, const std::vector<dso::TwoPartDate> &t,
            const std::vector<double> &y, int order, int &index,
            double &val) noexcept {
  constexpr const dso::DateTimeDifferenceType FD =
      dso::DateTimeDifferenceType::FractionalDays;

  // find suitable index in input dates
  if (index < 0) {
    auto it = std::lower_bound(t.begin(), t.end(), mjd);
    if (it == t.end())
      return 1;
    index = std::distance(t.begin(), it);
  }

  // make sure index is ok
  int window = (order + 1) / 2;
  if (index < window || index > (int)t.size() - window)
    return 1;

  // perform lagrangian interpolation, spanning indexes
  // [index-window, index+window)
  val = 0e0;
  for (int i = index - window; i < index + window; i++) {
    const double l = y[i];
    double pval = 1e0;
    for (int j = index - window; j < i; j++) {
      pval *= mjd.diff<FD>(t[j]) / t[i].diff<FD>(t[j]);
    }
    for (int j = i + 1; j < index + window; j++) {
      pval *= mjd.diff<FD>(t[j]) / t[i].diff<FD>(t[j]);
    }
    val += l * pval;
  }

  return 0;
}
} // unnamed namespace

int dso::EopLookUpTable::lagrange_interpolation(const dso::TwoPartDate &mjd,
                                                dso::EopRecord &eopr,
                                                int order) const noexcept {
  int status = 0;
  int index = -1;

  // xp (pole)
  status += lag_int(mjd, t, xp, order, index, eopr.xp);
  // yp (pole)
  status += lag_int(mjd, t, yp, order, index, eopr.yp);
  // Dut1
  status += lag_int(mjd, t, dut1, order, index, eopr.dut);
  // dX
  status += lag_int(mjd, t, dX, order, index, eopr.dx);
  // dY
  status += lag_int(mjd, t, dY, order, index, eopr.dy);
  // LOD
  status += lag_int(mjd, t, lod, order, index, eopr.lod);

  // remember to assign date to the filled-in instance
  eopr.mjd = mjd;

  return status;
}
