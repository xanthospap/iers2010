#include "cspline.hpp"
#include "datetime/dtfund.hpp"
#include "fundarg.hpp"
#include "hardisp.hpp"
#include "iau.hpp"
#include <algorithm>
#include <numeric>
#ifdef DEBUG
#include <cassert>
#endif

namespace {
struct Constituent {
  dso::DoodsonNumber d; ///< Doodson Number (no adding of +-5)
  double Hf;            ///< astronomical amplitude H_f [m]
};
// M2  S2  N2  K2  K1  O1  P1  Q1  MF  MM SSA
// Table 6.7 from IERS 2010
constexpr std::array<Constituent, dso::BlqSiteInfo::NTIN> BlqConstituents = {
    {/* M_2  */ {{2, 0, 0, 0, 0, 0}, .63192e0},
     /* S_2  */ {{2, 2, -2, 0, 0, 0}, .29400e0},
     /* N_2  */ {{2, -1, 0, 1, 0, 0}, .12099e0},
     /* K_2  */ {{2, 2, 0, 0, 0, 0}, .07996e0},
     /* K_1  */ {{1, 1, 0, 0, 0, 0}, .36878e0},
     /* O_1  */ {{1, -1, 0, 0, 0, 0}, -.26221e0},
     /* P_1  */ {{1, 1, -2, 0, 0, 0}, -.12203e0},
     /* Q_1  */ {{1, -2, 0, 1, 0, 0}, -.05020e0},
     /* M_f  */ {{0, 2, 0, 0, 0, 0}, -.06663e0},
     /* M_m  */ {{0, 1, 0, -1, 0, 0}, -.03518e0},
     /* S_sa */ {{0, 0, 2, 0, 0, 0}, -.03100e0}}};
constexpr const int NALL = 342;
constexpr std::array<Constituent, NALL> AllConstituents = {
    {{{2, 0, 0, 0, 0, 0}, +6.32208e-01},
     {{2, 2, -2, 0, 0, 0}, +2.94107e-01},
     {{2, -1, 0, 1, 0, 0}, +1.21046e-01},
     {{2, 2, 0, 0, 0, 0}, +7.99150e-02},
     {{2, 2, 0, 0, 1, 0}, +2.38180e-02},
     {{2, 0, 0, 0, -1, 0}, -2.35890e-02},
     {{2, -1, 2, -1, 0, 0}, +2.29940e-02},
     {{2, -2, 2, 0, 0, 0}, +1.93330e-02},
     {{2, 1, 0, -1, 0, 0}, -1.78710e-02},
     {{2, 2, -3, 0, 0, 1}, +1.71920e-02},
     {{2, -2, 0, 2, 0, 0}, +1.60180e-02},
     {{2, -3, 2, 1, 0, 0}, +4.67100e-03},
     {{2, 1, -2, 1, 0, 0}, -4.66200e-03},
     {{2, -1, 0, 1, -1, 0}, -4.51900e-03},
     {{2, 3, 0, -1, 0, 0}, +4.47000e-03},
     {{2, 1, 0, 1, 0, 0}, +4.46700e-03},
     {{2, 2, 0, 0, 2, 0}, +2.58900e-03},
     {{2, 2, -1, 0, 0, -1}, -2.45500e-03},
     {{2, 0, -1, 0, 0, 1}, -2.17200e-03},
     {{2, 1, 0, 1, 1, 0}, +1.97200e-03},
     {{2, 3, 0, -1, 1, 0}, +1.94700e-03},
     {{2, 0, 1, 0, 0, -1}, +1.91400e-03},
     {{2, 0, -2, 2, 0, 0}, -1.89800e-03},
     {{2, -3, 0, 3, 0, 0}, +1.80200e-03},
     {{2, -2, 3, 0, 0, -1}, +1.30400e-03},
     {{2, 4, 0, 0, 0, 0}, +1.17000e-03},
     {{2, -1, 1, 1, 0, -1}, +1.13000e-03},
     {{2, -1, 3, -1, 0, -1}, +1.06100e-03},
     {{2, 2, 0, 0, -1, 0}, -1.02200e-03},
     {{2, -1, -1, 1, 0, 1}, -1.01700e-03},
     {{2, 4, 0, 0, 1, 0}, +1.01400e-03},
     {{2, -3, 4, -1, 0, 0}, +9.01000e-04},
     {{2, -1, 2, -1, -1, 0}, -8.57000e-04},
     {{2, 3, -2, 1, 0, 0}, +8.55000e-04},
     {{2, 1, 2, -1, 0, 0}, +8.55000e-04},
     {{2, -4, 2, 2, 0, 0}, +7.72000e-04},
     {{2, 4, -2, 0, 0, 0}, +7.41000e-04},
     {{2, 0, 2, 0, 0, 0}, +7.41000e-04},
     {{2, -2, 2, 0, -1, 0}, -7.21000e-04},
     {{2, 2, -4, 0, 0, 2}, +6.98000e-04},
     {{2, 2, -2, 0, -1, 0}, +6.58000e-04},
     {{2, 1, 0, -1, -1, 0}, +6.54000e-04},
     {{2, -1, 1, 0, 0, 0}, -6.53000e-04},
     {{2, 2, -1, 0, 0, 1}, +6.33000e-04},
     {{2, 2, 1, 0, 0, -1}, +6.26000e-04},
     {{2, -2, 0, 2, -1, 0}, -5.98000e-04},
     {{2, -2, 4, -2, 0, 0}, +5.90000e-04},
     {{2, 2, 2, 0, 0, 0}, +5.44000e-04},
     {{2, -4, 4, 0, 0, 0}, +4.79000e-04},
     {{2, -1, 0, -1, -2, 0}, -4.64000e-04},
     {{2, 1, 2, -1, 1, 0}, +4.13000e-04},
     {{2, -1, -2, 3, 0, 0}, -3.90000e-04},
     {{2, 3, -2, 1, 1, 0}, +3.73000e-04},
     {{2, 4, 0, -2, 0, 0}, +3.66000e-04},
     {{2, 0, 0, 2, 0, 0}, +3.66000e-04},
     {{2, 0, 2, -2, 0, 0}, -3.60000e-04},
     {{2, 0, 2, 0, 1, 0}, -3.55000e-04},
     {{2, -3, 3, 1, 0, -1}, +3.54000e-04},
     {{2, 0, 0, 0, -2, 0}, +3.29000e-04},
     {{2, 4, 0, 0, 2, 0}, +3.28000e-04},
     {{2, 4, -2, 0, 1, 0}, +3.19000e-04},
     {{2, 0, 0, 0, 0, 2}, +3.02000e-04},
     {{2, 1, 0, 1, 2, 0}, +2.79000e-04},
     {{2, 0, -2, 0, -2, 0}, -2.74000e-04},
     {{2, -2, 1, 0, 0, 1}, -2.72000e-04},
     {{2, -2, 1, 2, 0, -1}, +2.48000e-04},
     {{2, -1, 1, -1, 0, 1}, -2.25000e-04},
     {{2, 5, 0, -1, 0, 0}, +2.24000e-04},
     {{2, 1, -3, 1, 0, 1}, -2.23000e-04},
     {{2, -2, -1, 2, 0, 1}, -2.16000e-04},
     {{2, 3, 0, -1, 2, 0}, +2.11000e-04},
     {{2, 1, -2, 1, -1, 0}, +2.09000e-04},
     {{2, 5, 0, -1, 1, 0}, +1.94000e-04},
     {{2, -4, 0, 4, 0, 0}, +1.85000e-04},
     {{2, -3, 2, 1, -1, 0}, -1.74000e-04},
     {{2, -2, 1, 1, 0, 0}, -1.71000e-04},
     {{2, 4, 0, -2, 1, 0}, +1.59000e-04},
     {{2, 0, 0, 2, 1, 0}, +1.31000e-04},
     {{2, -5, 4, 1, 0, 0}, +1.27000e-04},
     {{2, 0, 2, 0, 2, 0}, +1.20000e-04},
     {{2, -1, 2, 1, 0, 0}, +1.18000e-04},
     {{2, 5, -2, -1, 0, 0}, +1.17000e-04},
     {{2, 1, -1, 0, 0, 0}, +1.08000e-04},
     {{2, 2, -2, 0, 0, 2}, +1.07000e-04},
     {{2, -5, 2, 3, 0, 0}, +1.05000e-04},
     {{2, -1, -2, 1, -2, 0}, -1.02000e-04},
     {{2, -3, 5, -1, 0, -1}, +1.02000e-04},
     {{2, -1, 0, 0, 0, 1}, +9.90000e-05},
     {{2, -2, 0, 0, -2, 0}, -9.60000e-05},
     {{2, 0, -1, 1, 0, 0}, +9.50000e-05},
     {{2, -3, 1, 1, 0, 1}, -8.90000e-05},
     {{2, 3, 0, -1, -1, 0}, -8.50000e-05},
     {{2, 1, 0, 1, -1, 0}, -8.40000e-05},
     {{2, -1, 2, 1, 1, 0}, -8.10000e-05},
     {{2, 0, -3, 2, 0, 1}, -7.70000e-05},
     {{2, 1, -1, -1, 0, 1}, -7.20000e-05},
     {{2, -3, 0, 3, -1, 0}, -6.70000e-05},
     {{2, 0, -2, 2, -1, 0}, +6.60000e-05},
     {{2, -4, 3, 2, 0, -1}, +6.40000e-05},
     {{2, -1, 0, 1, -2, 0}, +6.30000e-05},
     {{2, 5, 0, -1, 2, 0}, +6.30000e-05},
     {{2, -4, 5, 0, 0, -1}, +6.30000e-05},
     {{2, -2, 4, 0, 0, -2}, +6.20000e-05},
     {{2, -1, 0, 1, 0, 2}, +6.20000e-05},
     {{2, -2, -2, 4, 0, 0}, -6.00000e-05},
     {{2, 3, -2, -1, -1, 0}, +5.60000e-05},
     {{2, -2, 5, -2, 0, -1}, +5.30000e-05},
     {{2, 0, -1, 0, -1, 1}, +5.10000e-05},
     {{2, 5, -2, -1, 1, 0}, +5.00000e-05},
     {{1, 1, 0, 0, 0, 0}, +3.68645e-01},
     {{1, -1, 0, 0, 0, 0}, -2.62232e-01},
     {{1, 1, -2, 0, 0, 0}, -1.21995e-01},
     {{1, -2, 0, 1, 0, 0}, -5.02080e-02},
     {{1, 1, 0, 0, 1, 0}, +5.00310e-02},
     {{1, -1, 0, 0, -1, 0}, -4.94700e-02},
     {{1, 2, 0, -1, 0, 0}, +2.06200e-02},
     {{1, 0, 0, 1, 0, 0}, +2.06130e-02},
     {{1, 3, 0, 0, 0, 0}, +1.12790e-02},
     {{1, -2, 2, -1, 0, 0}, -9.53000e-03},
     {{1, -2, 0, 1, -1, 0}, -9.46900e-03},
     {{1, -3, 2, 0, 0, 0}, -8.01200e-03},
     {{1, 0, 0, -1, 0, 0}, +7.41400e-03},
     {{1, 1, 0, 0, -1, 0}, -7.30000e-03},
     {{1, 3, 0, 0, 1, 0}, +7.22700e-03},
     {{1, 1, -3, 0, 0, 1}, -7.13100e-03},
     {{1, -3, 0, 2, 0, 0}, -6.64400e-03},
     {{1, 1, 2, 0, 0, 0}, +5.24900e-03},
     {{1, 0, 0, 1, 1, 0}, +4.13700e-03},
     {{1, 2, 0, -1, 1, 0}, +4.08700e-03},
     {{1, 0, 2, -1, 0, 0}, +3.94400e-03},
     {{1, 2, -2, 1, 0, 0}, +3.94300e-03},
     {{1, 3, -2, 0, 0, 0}, +3.42000e-03},
     {{1, -1, 2, 0, 0, 0}, +3.41800e-03},
     {{1, 1, 1, 0, 0, -1}, +2.88500e-03},
     {{1, 1, -1, 0, 0, 1}, +2.88400e-03},
     {{1, 4, 0, -1, 0, 0}, +2.16000e-03},
     {{1, -4, 2, 1, 0, 0}, -1.93600e-03},
     {{1, 0, -2, 1, 0, 0}, +1.93400e-03},
     {{1, -2, 2, -1, -1, 0}, -1.79800e-03},
     {{1, 3, 0, -2, 0, 0}, +1.69000e-03},
     {{1, -1, 0, 2, 0, 0}, +1.68900e-03},
     {{1, -1, 0, 0, -2, 0}, +1.51600e-03},
     {{1, 3, 0, 0, 2, 0}, +1.51400e-03},
     {{1, -3, 2, 0, -1, 0}, -1.51100e-03},
     {{1, 4, 0, -1, 1, 0}, +1.38300e-03},
     {{1, 0, 0, -1, -1, 0}, +1.37200e-03},
     {{1, 1, -2, 0, -1, 0}, +1.37100e-03},
     {{1, -3, 0, 2, -1, 0}, -1.25300e-03},
     {{1, 1, 0, 0, 2, 0}, -1.07500e-03},
     {{1, 1, -1, 0, 0, -1}, +1.02000e-03},
     {{1, -1, -1, 0, 0, 1}, +9.01000e-04},
     {{1, 0, 2, -1, 1, 0}, +8.65000e-04},
     {{1, -1, 1, 0, 0, -1}, -7.94000e-04},
     {{1, -1, -2, 2, 0, 0}, +7.88000e-04},
     {{1, 2, -2, 1, 1, 0}, +7.82000e-04},
     {{1, -4, 0, 3, 0, 0}, -7.47000e-04},
     {{1, -1, 2, 0, 1, 0}, -7.45000e-04},
     {{1, 3, -2, 0, 1, 0}, +6.70000e-04},
     {{1, 2, 0, -1, -1, 0}, -6.03000e-04},
     {{1, 0, 0, 1, -1, 0}, -5.97000e-04},
     {{1, -2, 2, 1, 0, 0}, +5.42000e-04},
     {{1, 4, -2, -1, 0, 0}, +5.42000e-04},
     {{1, -3, 3, 0, 0, -1}, -5.41000e-04},
     {{1, -2, 1, 1, 0, -1}, -4.69000e-04},
     {{1, -2, 3, -1, 0, -1}, -4.40000e-04},
     {{1, 0, -2, 1, -1, 0}, +4.38000e-04},
     {{1, -2, -1, 1, 0, 1}, +4.22000e-04},
     {{1, 4, -2, 1, 0, 0}, +4.10000e-04},
     {{1, -4, 4, -1, 0, 0}, -3.74000e-04},
     {{1, -4, 2, 1, -1, 0}, -3.65000e-04},
     {{1, 5, -2, 0, 0, 0}, +3.45000e-04},
     {{1, 3, 0, -2, 1, 0}, +3.35000e-04},
     {{1, -5, 2, 2, 0, 0}, -3.21000e-04},
     {{1, 2, 0, 1, 0, 0}, -3.19000e-04},
     {{1, 1, 3, 0, 0, -1}, +3.07000e-04},
     {{1, -2, 0, 1, -2, 0}, +2.91000e-04},
     {{1, 4, 0, -1, 2, 0}, +2.90000e-04},
     {{1, 1, -4, 0, 0, 2}, -2.89000e-04},
     {{1, 5, 0, -2, 0, 0}, +2.86000e-04},
     {{1, -1, 0, 2, 1, 0}, +2.75000e-04},
     {{1, -2, 1, 0, 0, 0}, +2.71000e-04},
     {{1, 4, -2, 1, 1, 0}, +2.63000e-04},
     {{1, -3, 4, -2, 0, 0}, -2.45000e-04},
     {{1, -1, 3, 0, 0, -1}, +2.25000e-04},
     {{1, 3, -3, 0, 0, 1}, +2.25000e-04},
     {{1, 5, -2, 0, 1, 0}, +2.21000e-04},
     {{1, 1, 2, 0, 1, 0}, -2.02000e-04},
     {{1, 2, 0, 1, 1, 0}, -2.00000e-04},
     {{1, -5, 4, 0, 0, 0}, -1.99000e-04},
     {{1, -2, 0, -1, -2, 0}, +1.92000e-04},
     {{1, 5, 0, -2, 1, 0}, +1.83000e-04},
     {{1, 1, 2, -2, 0, 0}, +1.83000e-04},
     {{1, 1, -2, 2, 0, 0}, +1.83000e-04},
     {{1, -2, 2, 1, 1, 0}, -1.70000e-04},
     {{1, 0, 3, -1, 0, -1}, +1.69000e-04},
     {{1, 2, -3, 1, 0, 1}, +1.68000e-04},
     {{1, -2, -2, 3, 0, 0}, +1.62000e-04},
     {{1, -1, 2, -2, 0, 0}, +1.49000e-04},
     {{1, -4, 3, 1, 0, -1}, -1.47000e-04},
     {{1, -4, 0, 3, -1, 0}, -1.41000e-04},
     {{1, -1, -2, 2, -1, 0}, +1.38000e-04},
     {{1, -2, 0, 3, 0, 0}, +1.36000e-04},
     {{1, 4, 0, -3, 0, 0}, +1.36000e-04},
     {{1, 0, 1, 1, 0, -1}, +1.27000e-04},
     {{1, 2, -1, -1, 0, 1}, +1.27000e-04},
     {{1, 2, -2, 1, -1, 0}, -1.26000e-04},
     {{1, 0, 0, -1, -2, 0}, -1.21000e-04},
     {{1, 2, 0, 1, 2, 0}, -1.21000e-04},
     {{1, 2, -2, -1, -1, 0}, +1.17000e-04},
     {{1, 0, 0, 1, 2, 0}, -1.16000e-04},
     {{1, 0, 1, 0, 0, 0}, -1.14000e-04},
     {{1, 2, -1, 0, 0, 0}, -1.14000e-04},
     {{1, 0, 2, -1, -1, 0}, -1.14000e-04},
     {{1, -1, -2, 0, -2, 0}, +1.14000e-04},
     {{1, -3, 1, 0, 0, 1}, +1.13000e-04},
     {{1, 3, -2, 0, -1, 0}, +1.09000e-04},
     {{1, -1, -1, 0, -1, 1}, +1.08000e-04},
     {{1, 4, -2, -1, 1, 0}, +1.06000e-04},
     {{1, 2, 1, -1, 0, -1}, -1.06000e-04},
     {{1, 0, -1, 1, 0, 1}, -1.06000e-04},
     {{1, -2, 4, -1, 0, 0}, +1.05000e-04},
     {{1, 4, -4, 1, 0, 0}, +1.04000e-04},
     {{1, -3, 1, 2, 0, -1}, -1.03000e-04},
     {{1, -3, 3, 0, -1, -1}, -1.00000e-04},
     {{1, 1, 2, 0, 2, 0}, -1.00000e-04},
     {{1, 1, -2, 0, -2, 0}, -1.00000e-04},
     {{1, 3, 0, 0, 3, 0}, +9.90000e-05},
     {{1, -1, 2, 0, -1, 0}, -9.80000e-05},
     {{1, -2, 1, -1, 0, 1}, +9.30000e-05},
     {{1, 0, -3, 1, 0, 1}, +9.30000e-05},
     {{1, -3, -1, 2, 0, 1}, +9.00000e-05},
     {{1, 2, 0, -1, 2, 0}, -8.80000e-05},
     {{1, 6, -2, -1, 0, 0}, +8.30000e-05},
     {{1, 2, 2, -1, 0, 0}, -8.30000e-05},
     {{1, -1, 1, 0, -1, -1}, -8.20000e-05},
     {{1, -2, 3, -1, -1, -1}, -8.10000e-05},
     {{1, -1, 0, 0, 0, 2}, -7.90000e-05},
     {{1, -5, 0, 4, 0, 0}, -7.70000e-05},
     {{1, 1, 0, 0, 0, -2}, -7.50000e-05},
     {{1, -2, 1, 1, -1, -1}, -7.50000e-05},
     {{1, 1, -1, 0, 1, 1}, -7.50000e-05},
     {{1, 1, 2, 0, 0, -2}, +7.10000e-05},
     {{1, -3, 1, 1, 0, 0}, +7.10000e-05},
     {{1, -4, 4, -1, -1, 0}, -7.10000e-05},
     {{1, 1, 0, -2, -1, 0}, +6.80000e-05},
     {{1, -2, -1, 1, -1, 1}, +6.80000e-05},
     {{1, -3, 2, 2, 0, 0}, +6.50000e-05},
     {{1, 5, -2, -2, 0, 0}, +6.50000e-05},
     {{1, 3, -4, 2, 0, 0}, +6.40000e-05},
     {{1, 1, -2, 0, 0, 2}, +6.40000e-05},
     {{1, -1, 4, -2, 0, 0}, +6.40000e-05},
     {{1, 2, 2, -1, 1, 0}, -6.40000e-05},
     {{1, -5, 2, 2, -1, 0}, -6.00000e-05},
     {{1, 1, -3, 0, -1, 1}, +5.60000e-05},
     {{1, 1, 1, 0, 1, -1}, +5.60000e-05},
     {{1, 6, -2, -1, 1, 0}, +5.30000e-05},
     {{1, -2, 2, -1, -2, 0}, +5.30000e-05},
     {{1, 4, -2, 1, 2, 0}, +5.30000e-05},
     {{1, -6, 4, 1, 0, 0}, -5.30000e-05},
     {{1, 5, -4, 0, 0, 0}, +5.30000e-05},
     {{1, -3, 4, 0, 0, 0}, +5.30000e-05},
     {{1, 1, 2, -2, 1, 0}, +5.20000e-05},
     {{1, -2, 1, 0, -1, 0}, +5.00000e-05},
     {{0, 2, 0, 0, 0, 0}, -6.66070e-02},
     {{0, 1, 0, -1, 0, 0}, -3.51840e-02},
     {{0, 0, 2, 0, 0, 0}, -3.09880e-02},
     {{0, 0, 0, 0, 1, 0}, +2.79290e-02},
     {{0, 2, 0, 0, 1, 0}, -2.76160e-02},
     {{0, 3, 0, -1, 0, 0}, -1.27530e-02},
     {{0, 1, -2, 1, 0, 0}, -6.72800e-03},
     {{0, 2, -2, 0, 0, 0}, -5.83700e-03},
     {{0, 3, 0, -1, 1, 0}, -5.28600e-03},
     {{0, 0, 1, 0, 0, -1}, -4.92100e-03},
     {{0, 2, 0, -2, 0, 0}, -2.88400e-03},
     {{0, 2, 0, 0, 2, 0}, -2.58300e-03},
     {{0, 3, -2, 1, 0, 0}, -2.42200e-03},
     {{0, 1, 0, -1, -1, 0}, +2.31000e-03},
     {{0, 1, 0, -1, 1, 0}, +2.28300e-03},
     {{0, 4, -2, 0, 0, 0}, -2.03700e-03},
     {{0, 1, 0, 1, 0, 0}, +1.88300e-03},
     {{0, 0, 3, 0, 0, -1}, -1.81100e-03},
     {{0, 4, 0, -2, 0, 0}, -1.68700e-03},
     {{0, 3, -2, 1, 1, 0}, -1.00400e-03},
     {{0, 3, -2, -1, 0, 0}, -9.25000e-04},
     {{0, 4, -2, 0, 1, 0}, -8.44000e-04},
     {{0, 0, 2, 0, 1, 0}, +7.66000e-04},
     {{0, 1, 0, 1, 1, 0}, +7.66000e-04},
     {{0, 4, 0, -2, 1, 0}, -7.00000e-04},
     {{0, 3, 0, -1, 2, 0}, -4.95000e-04},
     {{0, 5, -2, -1, 0, 0}, -4.92000e-04},
     {{0, 1, 2, -1, 0, 0}, +4.91000e-04},
     {{0, 1, -2, 1, -1, 0}, +4.83000e-04},
     {{0, 1, -2, 1, 1, 0}, +4.37000e-04},
     {{0, 2, -2, 0, -1, 0}, -4.16000e-04},
     {{0, 2, -3, 0, 0, 1}, -3.84000e-04},
     {{0, 2, -2, 0, 1, 0}, +3.74000e-04},
     {{0, 0, 2, -2, 0, 0}, -3.12000e-04},
     {{0, 1, -3, 1, 0, 1}, -2.88000e-04},
     {{0, 0, 0, 0, 2, 0}, -2.73000e-04},
     {{0, 0, 1, 0, 0, 1}, +2.59000e-04},
     {{0, 1, 2, -1, 1, 0}, +2.45000e-04},
     {{0, 3, 0, -3, 0, 0}, -2.32000e-04},
     {{0, 2, 1, 0, 0, -1}, +2.29000e-04},
     {{0, 1, -1, -1, 0, 1}, -2.16000e-04},
     {{0, 1, 0, 1, 2, 0}, +2.06000e-04},
     {{0, 5, -2, -1, 1, 0}, -2.04000e-04},
     {{0, 2, -1, 0, 0, 1}, -2.02000e-04},
     {{0, 2, 2, -2, 0, 0}, +2.00000e-04},
     {{0, 1, -1, 0, 0, 0}, +1.95000e-04},
     {{0, 5, 0, -3, 0, 0}, -1.90000e-04},
     {{0, 2, 0, -2, 1, 0}, +1.87000e-04},
     {{0, 1, 1, -1, 0, -1}, +1.80000e-04},
     {{0, 3, -4, 1, 0, 0}, -1.79000e-04},
     {{0, 0, 2, 0, 2, 0}, +1.70000e-04},
     {{0, 2, 0, -2, -1, 0}, +1.53000e-04},
     {{0, 4, -3, 0, 0, 1}, -1.37000e-04},
     {{0, 3, -1, -1, 0, 1}, -1.19000e-04},
     {{0, 0, 2, 0, 0, -2}, -1.19000e-04},
     {{0, 3, -3, 1, 0, 1}, -1.12000e-04},
     {{0, 2, -4, 2, 0, 0}, -1.10000e-04},
     {{0, 4, -2, -2, 0, 0}, -1.10000e-04},
     {{0, 3, 1, -1, 0, -1}, +1.07000e-04},
     {{0, 5, -4, 1, 0, 0}, -9.50000e-05},
     {{0, 3, -2, -1, -1, 0}, -9.50000e-05},
     {{0, 3, -2, 1, 2, 0}, -9.10000e-05},
     {{0, 4, -4, 0, 0, 0}, -9.00000e-05},
     {{0, 6, -2, -2, 0, 0}, -8.10000e-05},
     {{0, 5, 0, -3, 1, 0}, -7.90000e-05},
     {{0, 4, -2, 0, 2, 0}, -7.90000e-05},
     {{0, 2, 2, -2, 1, 0}, +7.70000e-05},
     {{0, 0, 4, 0, 0, -2}, -7.30000e-05},
     {{0, 3, -1, 0, 0, 0}, +6.90000e-05},
     {{0, 3, -3, -1, 0, 1}, -6.70000e-05},
     {{0, 4, 0, -2, 2, 0}, -6.60000e-05},
     {{0, 1, -2, -1, -1, 0}, +6.50000e-05},
     {{0, 2, -1, 0, 0, -1}, +6.40000e-05},
     {{0, 4, -4, 2, 0, 0}, -6.20000e-05},
     {{0, 2, 1, 0, 1, -1}, +6.00000e-05},
     {{0, 3, -2, -1, 1, 0}, +5.90000e-05},
     {{0, 4, -3, 0, 1, 1}, -5.60000e-05},
     {{0, 2, 0, 0, 3, 0}, +5.50000e-05},
     {{0, 6, -4, 0, 0, 0}, -5.10000e-05}}};

} // unnamed namespace

/// @brief Compute (Doodson) angles beta and beta_freq for given epoch.
/// These angles enable the computation of frequency and phase for any
/// tidal constituent (at given epoch)
int dso::Hardisp::operator()(const dso::TwoPartDate &tt_mjd) noexcept {
  // compute angles:
  // A1. Fundamental Arguments (temporary)
  double fundarg[5];
  iers2010::fundarg(tt_mjd, fundarg);
  // A2. GMST [rad]
  const double gmst = /*iers2010::sofa::gmst06(ut1_jd._big, ut1_jd._small,
                                             tt_jd._big, tt_jd._small);*/
      gmst_utc(tt_mjd.tt2utc());
  // A3. Doodson variables (stored in beta)
  dso::fundarg2doodson(fundarg, gmst, beta);
  // TODO correct first Doodson variable to match TDFRPH
  // D(1) = 360.0D0*DAYFR - F4
  // beta[0] = ut1_mjd._small * dso::D2PI - fundarg[3];

  dso::doodson_frequency_args(tt_mjd, beta_freq);
  // TODO why the fuck does this not work ??
  // B1. Derivatives of Fundamental Arguments in [rad/century]
  // iers2010::fundarg_derivs(tt_mjd.jcenturies_sinceJ2000(), fundarg);
  // B2. Derivative of GMST angle [rad/century]
  // const double dgmst = iers2010::sofa::dgmst06dt(tt_jd._big, tt_jd._small);
  // B3. Variables for frequency (stored in in beta_freq)
  // dso::fundarg_derivs2doodson_derivs(fundarg, dgmst, beta_freq);
  // B4. [rad/century] to [cpd]
  // for (int i=0; i<6; i++) beta_freq[i] /= (dso::days_in_julian_cent /
  // dso::D2PI);
  // tdfrph(fundarg, beta_freq);
  // for (int i=0; i<6; i++) printf("\tdd(%d)%+.6f\n", i+1, beta[i]);

  //#ifdef DEBUG
  //  const double tp = tt_mjd.jcenturies_sinceJ2000();
  //  printf("Time argument in TDFRPH is: %.12f\n", tp);
  //#endif
  return 0;
}

void dso::Hardisp::hardisp(const dso::BlqSiteInfo &blq, double &dr, double &dw,
                           double &ds) noexcept {
  // call admint function and store results in instance's struct
  admint(blq);
  // normalized frequencies; cannot be over NALL
  double scr[NALL * 3 * 3];
  double *__restrict__ scrr = scr;
  double *__restrict__ scrw = scr + 3 * NALL;
  double *__restrict__ scrs = scr + 6 * NALL;
  // normalize frequencies and c
  for (int i = 0; i < admnt.num_constituents(); i++) {
    const int j = 3 * i;
    const double wf = admnt[i].freq; // * dso::D2PI / dso::sec_per_day;
    // radial component
    const double rc = admnt[i].ampl_r * std::cos(admnt[i].ph_r);
    const double rs = -admnt[i].ampl_r * std::sin(admnt[i].ph_r);
    // printf("Constituent %2d wf=%.6f rc=%.6f rs=%.6f\n", i, wf, rc, rs);
    scrr[j + 0] = rc;
    // printf("scrr = %.6f\n", scrr[i]);
    scrr[j + 1] = rc * std::cos(wf) - rs * std::sin(wf);
    scrr[j + 2] = 2e0 * std::cos(wf);
    // east-west component
    const double wc = admnt[i].ampl_w * std::cos(admnt[i].ph_w);
    const double ws = -admnt[i].ampl_w * std::sin(admnt[i].ph_w);
    scrw[j + 0] = wc;
    scrw[j + 1] = wc * std::cos(wf) - ws * std::sin(wf);
    scrw[j + 2] = 2e0 * std::cos(wf);
    // north-south component
    const double sc = admnt[i].ampl_s * std::cos(admnt[i].ph_s);
    const double ss = -admnt[i].ampl_s * std::sin(admnt[i].ph_s);
    scrs[j + 0] = sc;
    scrs[j + 1] = sc * std::cos(wf) - ss * std::sin(wf);
    scrs[j + 2] = 2e0 * std::cos(wf);
  }
  double r = 0e0, w = 0e0, s = 0e0;
  for (int i = 0; i < admnt.num_constituents(); i++) {
    const int j = 3 * i;
    r += scrr[j];
    // printf("X + %12.6f\n", scrr[j]);
    double sc = scrr[j];
    scrr[j] = scrr[j + 2] * sc - scrr[j + 1];
    scrr[j + 1] = sc;
    w += scrw[j];
    sc = scrw[j];
    scrw[j] = scrw[j + 2] * sc - scrw[j + 1];
    scrw[j + 1] = sc;
    s += scrs[j];
    sc = scrs[j];
    scrs[j] = scrs[j + 2] * sc - scrs[j + 1];
    scrs[j + 1] = sc;
  }
  dr = r;
  dw = w;
  ds = s;
  return;
}

void dso::Hardisp::admint(const dso::BlqSiteInfo &blq) noexcept {
  constexpr const int N = dso::BlqSiteInfo::NTIN;
  // indexes of AllConstituents matching the ones in dso::BlqSiteInfo::Doodsons
  double pool[7 * N];
  double *__restrict__ rf = pool;
  double *__restrict__ aim = rf + N;     ///< radial, west, south
  double *__restrict__ rl = aim + 3 * N; ///< radial, west, south
  // char buf[32];

  admnt.clear();
  admnt.reserve(NALL);

  // match given Doodson numbers to the ones in AllConstituents and compute
  // in radial, west, south components
  // Here, we have collected the 11 fundamental constituents of BLQ, so as not
  // to search through AllConstituents
  int i = 0;
  for (const auto &c : BlqConstituents) {
    // RL(K) = AMPIN(LL)*COS(DTR*PHIN(LL))/ABS(TAMP(KK))
    rl[i] = blq.amplitudes[i] * std::cos(-blq.phases[i]) / std::abs(c.Hf);
    // printf("\tRL  = %.9f\n", rl[i]);
    rl[i + N] =
        blq.amplitudes[i + N] * std::cos(-blq.phases[i + N]) / std::abs(c.Hf);
    rl[i + 2 * N] = blq.amplitudes[i + 2 * N] *
                    std::cos(-blq.phases[i + 2 * N]) / std::abs(c.Hf);
    // AIM(K)= AMPIN(LL)*SIN(DTR*PHIN(LL))/ABS(TAMP(KK))
    aim[i] = blq.amplitudes[i] * std::sin(-blq.phases[i]) / std::abs(c.Hf);
    // printf("\tAIM = %.9f\n", aim[i]);
    aim[i + N] =
        blq.amplitudes[i + N] * std::sin(-blq.phases[i + N]) / std::abs(c.Hf);
    aim[i + 2 * N] = blq.amplitudes[i + 2 * N] *
                     std::sin(-blq.phases[i + 2 * N]) / std::abs(c.Hf);
    // frequency [rad/century] to [cpd]
    // rf[i] = c.d.frequency(beta_freq) / dso::days_in_julian_cent / dso::D2PI;
    rf[i] = c.d.frequency(beta_freq) / dso::D2PI;
    // printf("\tRF  = %.9f\n", rf[i]);
    //  augment counter
    ++i;
  }

  // sort constituents based on RF array (aka frequency)
  {
    // use k array to hold sorted indexes, based on rf
    std::array<int, N> k;
    std::iota(k.begin(), k.end(), 0);
    std::sort(k.begin(), k.end(),
              [rf](int i1, int i2) { return *(rf + i1) < *(rf + i2); });
    // based on k, sort/re-order arrays rl and aim, per component
    double tmp[N];
    for (int j = 0; j < 3; j++) {
      const int offset = j * N;
      std::memcpy(tmp, rl + offset, sizeof(double) * N);
      for (int l = 0; l < N; l++)
        *(rl + offset + l) = tmp[k[l]];
    }
    for (int j = 0; j < 3; j++) {
      const int offset = j * N;
      std::memcpy(tmp, aim + offset, sizeof(double) * N);
      for (int l = 0; l < N; l++)
        *(aim + offset + l) = tmp[k[l]];
    }
    // last thing, reorder rf (no components, single array)
    std::memcpy(tmp, rf, sizeof(double) * N);
    for (int l = 0; l < N; l++)
      *(rf + l) = tmp[k[l]];
  }
  // for (int j = 0; j < N; j++) {
  //   printf("Constituent %2d: rf: %+.6f rl: %+.6f aim: %+.6f\n",
  //       j + 1, rf[j], rl[j], aim[j]);
  // }

  // indexes/sizes for spline interpolation based on frequencies
  int nlp = 0;
  int ndi = 0;
  int nsd = 0;
  for (i = 0; i < N; i++) {
    if (rf[i] < .5e0)
      ++nlp;
    else if (rf[i] < 1.5e0)
      ++ndi;
    else if (rf[i] < 2.5e0)
      ++nsd;
  }

  double scr[N];
  double d2[N * 6];
  // one component at a time, for radial, west, south
  for (i = 0; i < 3; i++) {
    double *__restrict__ caim = aim + i * N; ///< radial, west, south
    double *__restrict__ crl = rl + i * N;   ///< radial, west, south
    // setup arrays for spline interpolation
    if (nlp) {
      // IF(NLP.NE.0) CALL SPLINE(NLP,RF,RL,ZDR,SCR)
      cspline_deriv(rf, crl, nlp, d2 + 0 * N, 0, 0, scr);
      // IF(NLP.NE.0) CALL SPLINE(NLP,RF,AIM,ZDI,SCR)
      cspline_deriv(rf, caim, nlp, d2 + 1 * N, 0, 0, scr);
    }
    cspline_deriv(rf + nlp, crl + nlp, ndi, d2 + 2 * N, 0, 0, scr);
    cspline_deriv(rf + nlp, caim + nlp, ndi, d2 + 3 * N, 0, 0, scr);
    cspline_deriv(rf + nlp + ndi, crl + nlp + ndi, nsd, d2 + 4 * N, 0, 0, scr);
    cspline_deriv(rf + nlp + ndi, caim + nlp + ndi, nsd, d2 + 5 * N, 0, 0, scr);
    // printf("spline with n=%d starting from %d\n", nlp, 0);
    // printf("spline with n=%d starting from %d\n", ndi, nlp);
    // printf("spline with n=%d starting from %d\n", nsd, nlp+ndi);

    int index = 0;
    // Evaluate all harmonics using the interpolated admittance
    for (const auto &cnst : AllConstituents) {
      // IF(IDD(1,I).EQ.0.AND.NLP.EQ.0) GO TO 11
      if (!(!cnst.d(0) && !nlp)) {
        // CALL TDFRPH(IDD(1,I),F(J),P(J))
        const double phase = cnst.d.phase(beta) + (cnst.d(0) == 0) * dso::DPI +
                             (cnst.d(0) == 1) * (dso::DPI / 2e0);
        const double freq = cnst.d.frequency(beta_freq) / dso::D2PI;
        // printf("interpolation for const=%d at freq=%.6f, phase=%.6f\n",
        // (int)(&cnst-&AllConstituents[0]),freq,dso::rad2deg(phase));
        int num_pts = 0, ar_offset = 0, d2_offset = 0;
        switch (cnst.d(0)) {
        case 0:
          num_pts = nlp;
          d2_offset = ar_offset = 0;
          break;
        case 1:
          num_pts = ndi;
          ar_offset = nlp;
          d2_offset = 2;
          break;
        case 2:
          num_pts = nsd;
          ar_offset = nlp + ndi;
          d2_offset = 4;
        }
        int guess = -1;
        double re, am;
        cspline_interp(freq, guess, rf + ar_offset, crl + ar_offset, num_pts,
                       d2 + d2_offset * N, re);
        cspline_interp(freq, guess, rf + ar_offset, caim + ar_offset, num_pts,
                       d2 + (d2_offset + 1) * N, am);
        // printf("RE=%15.6f AM=%15.6f D(0)=%d num_pts=%d start=%d\n", re, am,
        // cnst.d(0), num_pts, ar_offset);
        const double amp = cnst.Hf * std::sqrt(re * re + am * am);
        const double ph = dso::anpm(std::atan2(am, re) + phase);
        admnt.add_component(freq, amp, ph, i, index);
        ++index;
        if (false)
          printf(
              "Constituent %d Amplitude: %+12.9f Phase: %+12.6f Freq: %+.9f\n",
              (int)(&cnst - &AllConstituents[0]), amp, dso::rad2deg(ph), freq);
      }
    }
  }

  return;
}
