#include "eop.hpp"

std::vector<double> dso::EopLookUpTable::dutToUt1Tai() const noexcept {
  std::vector<double> dUt1Tai (dut1);
  for (auto epoch = t.cbegin(); epoch != t.cend(); ++epoch) {
    // index
    const int i = std::distance(t.cbegin(), epoch);
    // need ΔAT for this UTC
    const dso::modified_julian_day mjd((int)(epoch->as_mjd()));
    const int dat = dso::dat(mjd);
    // construct [UT1-UTC] = [UT1-UTC] - ΔAT
    dUt1Tai[i] -= static_cast<double>(dat);
  }
  return dUt1Tai;
}
