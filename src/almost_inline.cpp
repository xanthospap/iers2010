#include "iers2010.hpp"

int iers2010::ortho_eop(double fmjd, double &dx, double &dy,
                        double &dut1) noexcept {
  return ortho_eop(split_fmjd(fmjd), dx, dy, dut1);
}
int iers2010::oeop::cnmtx(double fmjd, double *h) noexcept {
  return cnmtx(split_fmjd(fmjd), h);
}
