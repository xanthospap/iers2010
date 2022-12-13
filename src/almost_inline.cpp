#include "iers2010.hpp"

int iers2010::pmsdnut2(double fmjd, double &dx, double &dy) noexcept {
  return pmsdnut2(split_fmjd(fmjd), dx, dy);
}
int iers2010::utlibr(double fmjd, double &dut1, double &dlod) noexcept {
  return utlibr(split_fmjd(fmjd), dut1, dlod);
}
int iers2010::ortho_eop(double fmjd, double &dx, double &dy,
                        double &dut1) noexcept {
  return ortho_eop(split_fmjd(fmjd), dx, dy, dut1);
}
int iers2010::oeop::cnmtx(double fmjd, double *h) noexcept {
  return cnmtx(split_fmjd(fmjd), h);
}
