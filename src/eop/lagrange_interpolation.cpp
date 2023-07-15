#include "eop.hpp"
#include <algorithm>

namespace {
int lag_int(const dso::TwoPartDate &mjd, const std::vector<dso::TwoPartDate> &t,
            const std::vector<double> &y, int order, int &index,
            double &val) noexcept {

  /* Units of delta date; use fractional days */
  constexpr const dso::DateTimeDifferenceType FD =
      dso::DateTimeDifferenceType::FractionalDays;

  /* find suitable index in input dates if the passed in index is < 0 */
  if (index < 0) {
    /* it->t > mjd */
    auto it = std::upper_bound(t.begin(), t.end(), mjd);
    /* mjd is past the last date in the the t vector */
    if (it == t.end())
      return 1;
    index = std::distance(t.begin(), it);
  }

  /* make sure we have enough data points */
  int window = (order + 1) / 2;
  int left_start = index - window + 1;
  int right_stop = index + window;
  if (left_start < 0 || right_stop >= (int)t.size())
    return 1;

  /* perform lagrangian interpolation, spanning indexes 
   * [index-window+1, index+window]
   */
  val = 0e0;
  for (int i = left_start; i < right_stop + 1; i++) {
    const double l = y[i];
    double pval = 1e0;
    for (int j = left_start; j < i; j++) {
      pval *= mjd.diff<FD>(t[j]) / t[i].diff<FD>(t[j]);
    }
    for (int j = i + 1; j < right_stop + 1; j++) {
      pval *= mjd.diff<FD>(t[j]) / t[i].diff<FD>(t[j]);
    }
    val += l * pval;
  }

  return 0;
}
} /* unnamed namespace */

int dso::EopLookUpTable::lagrange_interpolation(const dso::TwoPartDate &mjd,
                                                dso::EopRecord &eopr,
                                                int order) const noexcept {
  int status = 0;

  /* first time in lag_int use a negative index to signal that we are 
   * searching for the right interval; next calls will use the index already 
   * obtained from the first call
   */
  int index = -1;

  /* xp (pole) */
  status += lag_int(mjd, tvec, xp, order, index, eopr.xp);
  /* yp (pole) */
  status += lag_int(mjd, tvec, yp, order, index, eopr.yp);
  /* Dut1 */
  status += lag_int(mjd, tvec, dut1, order, index, eopr.dut);
  /* dX */
  status += lag_int(mjd, tvec, dX, order, index, eopr.dx);
  /* dY */
  status += lag_int(mjd, tvec, dY, order, index, eopr.dy);
  /* LOD */
  status += lag_int(mjd, tvec, lod, order, index, eopr.lod);

  /* remember to assign date to the filled-in instance */
  eopr.mjd = mjd;

  return status;
}
