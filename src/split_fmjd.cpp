#include "iers2010.hpp"

dso::TwoPartDate iers2010::split_fmjd(double fmjd) noexcept {
  double ip;
  const double fp = std::modf(fmjd, &ip);
  return dso::TwoPartDate(ip,fp);
}
