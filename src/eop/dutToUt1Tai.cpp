#include "eop.hpp"

std::vector<double> dso::EopLookUpTable::dutToUt1Tai() const noexcept {
  /* copy of original vector (ΔUT1) */
  std::vector<double> dUt1Tai (dut1);

  for (auto epoch = tvec.cbegin(); epoch != tvec.cend(); ++epoch) {
    /* current index */
    const int i = std::distance(tvec.cbegin(), epoch);
    /* need ΔAT for this UTC */
    const dso::modified_julian_day imjd((long)(epoch->as_mjd()));
    const int dat = dso::dat(imjd);
    /* construct [UT1-UTC] = [UT1-UTC] - ΔAT */
    dUt1Tai[i] -= static_cast<double>(dat);
  }

  /* return vector */
  return dUt1Tai;
}
