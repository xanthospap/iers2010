#include "doodson.hpp"
#include "eop.hpp"
#include <array>
#include <cmath>

namespace {
struct Table51aData {
  /* Doodson numbers from IERS 2010, Table 5.1a */
  dso::DoodsonConstituent dn;
  /* cooresponding Coefficients of sin(argument) and cos(argument),
   * IERS 2010, Table 5.1a
   */
  double xs, xc, ys, yc;
};
constexpr const std::array<Table51aData, 11> Table = {{
    {{+2, -2, +0, +2, +0, +0},
     +5.000e-02,
     -3.000e-02,
     -3.000e-01,
     -6.000e-01}, /*      0.54*/
    {{+2, -2, +2, +0, +0, +0},
     +6.000e-02,
     -3.000e-02,
     -4.000e-01,
     -7.000e-01}, /*      0.54*/
    {{+2, -1, +0, +1, +0, +0},
     +3.500e-01,
     -2.000e-01,
     -2.400e+00,
     -4.200e+00}, /*      0.53*/
    {{+2, -1, +2, -1, +0, +0},
     +7.000e-02,
     -4.000e-02,
     -5.000e-01,
     -8.000e-01}, /*      0.53*/
    {{+2, +0, +0, +0, -1, +0},
     -7.000e-02,
     +4.000e-02,
     +5.000e-01,
     +8.000e-01}, /*      0.52*/
    {{+2, +0, +0, +0, +0, +0},
     +1.750e+00,
     -1.010e+00,
     -1.220e+01,
     -2.130e+01}, /*      0.52*/
    {{+2, +1, +0, -1, +0, +0},
     -5.000e-02,
     +3.000e-02,
     +3.000e-01,
     +6.000e-01}, /*      0.51*/
    {{+2, +2, -3, +0, +0, +1},
     +5.000e-02,
     -3.000e-02,
     -3.000e-01,
     -6.000e-01}, /*      0.50*/
    {{+2, +2, -2, +0, +0, +0},
     +7.600e-01,
     -4.400e-01,
     -5.500e+00,
     -9.500e+00}, /*      0.50*/
    {{+2, +2, +0, +0, +0, +0},
     +2.100e-01,
     -1.200e-01,
     -1.500e+00,
     -2.600e+00}, /*      0.50*/
    {{+2, +2, +0, +0, +1, +0},
     +6.000e-02,
     -4.000e-02,
     -4.000e-01,
     -8.000e-01}, /*      0.50*/
}};
} /* unnamed namespace */

int dso::utlod_libration(const double *const fargs, double gmst, double &dut1,
                     double &dlod) noexcept {
  /* get Doodson arguments from Delaunay arguments and GMST */
  double dargs[6];
  dso::delaunay2doodson(fargs, gmst, dargs);
  /* these will never change */
  const double *__restrict__ f = dargs;

  /* set corrections to zero */
  dlod = dut1 = 0e0;
  for (const auto &entry : Table) {
    const double arg = entry.dn.argument(f);
    const double sa = std::sin(arg);
    const double ca = std::cos(arg);
    dut1 += sa * entry.xs + ca * entry.xc; /* microseconds */
    dlod += sa * entry.ys + ca * entry.yc; /* microseconds per day */
  }

  return 0;
}
