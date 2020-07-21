#include "iers2010.hpp"
#include <cassert>
#ifdef USE_EXTERNAL_CONSTS
#include "gencon.hpp"
#endif

/// @details  The purpose of the subroutine is to compute the time dependent
///           part of second degree diurnal and semidiurnal tidal potential
///           from the dominant spectral lines in the Cartwright-Tayler-Edden
///           harmonic decomposition.
///           This function is a translation/wrapper for the fortran CNMTX
///           subroutine, found here :
///           http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]   dmjd  Modified Julian Date
/// @param[out]  h     Vector of length 12 with partials of the tidal variation
///                    with respect to the orthoweights (Note 1)
/// @return            An integer, always 0.
///
/// @note
///     -# The diurnal and semidiurnal orthoweights fit to the 8 constituents
///        are listed in Reference Ray et al.
///     -# Status: Canonical model
///
/// @version 17.03.2010
///
/// @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS
///       Technical Note No. 36, BKG (2010)
/// Ray,R. D., Steinberg, D. J., Chao, B. F., and Cartwright, D. E.,
///      "Diurnal and Semidiurnal Variations in the Earth's Rotation
///      Rate Induced by Ocean Tides", 1994, Science, 264, pp. 830-832
int iers2010::oeop::cnmtx(double dmjd, double *h) {
#ifdef USE_EXTERNAL_CONSTS
  constexpr double TWOPI(D2PI);
#else
  constexpr double TWOPI(6.283185307179586476925287e0);
#endif

  // Define the orthotide weight factors
  constexpr double sp[2][6] = {
      {0.0298e0, 0.1408e0, +0.0805e0, 0.6002e0, +0.3025e0, 0.1517e0},
      {0.0200e0, 0.0905e0, +0.0638e0, 0.3476e0, +0.1645e0, 0.0923e0}};

  constexpr double dt(2e0);
  constexpr int nmax(2);

  // tidal potential model for 71 diurnal and semidiurnal lines
  constexpr double d1960(37076.5e0);
  constexpr struct {
    int nj, mj;
    double hs, phase, freq;
    char numarg[8];
  } x[] = {// DATA (NJ(J),MJ(J),HS(J),PHASE(J),FREQ(J),NUMARG(J),J=1,15)
           {2, 1, -1.94e0, 9.0899831e0, 5.18688050e0, "117.655"},
           {2, 1, -1.25e0, 8.8234208e0, 5.38346657e0, "125.745"},
           {2, 1, -6.64e0, 12.1189598e0, 5.38439079e0, "125.755"},
           {2, 1, -1.51e0, 1.4425700e0, 5.41398343e0, "127.545"},
           {2, 1, -8.02e0, 4.7381090e0, 5.41490765e0, "127.555"},
           {2, 1, -9.47e0, 4.4715466e0, 5.61149372e0, "135.645"},
           {2, 1, -50.20e0, 7.7670857e0, 5.61241794e0, "135.655"},
           {2, 1, -1.80e0, -2.9093042e0, 5.64201057e0, "137.445"},
           {2, 1, -9.54e0, 0.3862349e0, 5.64293479e0, "137.455"},
           {2, 1, 1.52e0, -3.1758666e0, 5.83859664e0, "145.535"},
           {2, 1, -49.45e0, 0.1196725e0, 5.83952086e0, "145.545"},
           {2, 1, -262.21e0, 3.4152116e0, 5.84044508e0, "145.555"},
           {2, 1, 1.70e0, 12.8946194e0, 5.84433381e0, "145.755"},
           {2, 1, 3.43e0, 5.5137686e0, 5.87485066e0, "147.555"},
           {2, 1, 1.94e0, 6.4441883e0, 6.03795537e0, "153.655"},
           // DATA (NJ(J),MJ(J),HS(J),PHASE(J),FREQ(J),NUMARG(J),J=16,30)
           {2, 1, 1.37e0, -4.2322016e0, 6.06754801e0, "155.445"},
           {2, 1, 7.41e0, -0.9366625e0, 6.06847223e0, "155.455"},
           {2, 1, 20.62e0, 8.5427453e0, 6.07236095e0, "155.655"},
           {2, 1, 4.14e0, 11.8382843e0, 6.07328517e0, "155.665"},
           {2, 1, 3.94e0, 1.1618945e0, 6.10287781e0, "157.455"},
           {2, 1, -7.14e0, 5.9693878e0, 6.24878055e0, "162.556"},
           {2, 1, 1.37e0, -1.2032249e0, 6.26505830e0, "163.545"},
           {2, 1, -122.03e0, 2.0923141e0, 6.26598252e0, "163.555"},
           {2, 1, 1.02e0, -1.7847596e0, 6.28318449e0, "164.554"},
           {2, 1, 2.89e0, 8.0679449e0, 6.28318613e0, "164.556"},
           {2, 1, -7.30e0, 0.8953321e0, 6.29946388e0, "165.545"},
           {2, 1, 368.78e0, 4.1908712e0, 6.30038810e0, "165.555"},
           {2, 1, 50.01e0, 7.4864102e0, 6.30131232e0, "165.565"},
           {2, 1, -1.08e0, 10.7819493e0, 6.30223654e0, "165.575"},
           {2, 1, 2.93e0, 0.3137975e0, 6.31759007e0, "166.554"},
           // DATA (NJ(J),MJ(J),HS(J),PHASE(J),FREQ(J),NUMARG(J),J=31,45)
           {2, 1, 5.25e0, 6.2894282e0, 6.33479368e0, "167.555"},
           {2, 1, 3.95e0, 7.2198478e0, 6.49789839e0, "173.655"},
           {2, 1, 20.62e0, -0.1610030e0, 6.52841524e0, "175.455"},
           {2, 1, 4.09e0, 3.1345361e0, 6.52933946e0, "175.465"},
           {2, 1, 3.42e0, 2.8679737e0, 6.72592553e0, "183.555"},
           {2, 1, 1.69e0, -4.5128771e0, 6.75644239e0, "185.355"},
           {2, 1, 11.29e0, 4.9665307e0, 6.76033111e0, "185.555"},
           {2, 1, 7.23e0, 8.2620698e0, 6.76125533e0, "185.565"},
           {2, 1, 1.51e0, 11.5576089e0, 6.76217955e0, "185.575"},
           {2, 1, 2.16e0, 0.6146566e0, 6.98835826e0, "195.455"},
           {2, 1, 1.38e0, 3.9101957e0, 6.98928248e0, "195.465"},
           {2, 2, 1.80e0, 20.6617051e0, 11.45675174e0, "225.855"},
           {2, 2, 4.67e0, 13.2808543e0, 11.48726860e0, "227.655"},
           {2, 2, 16.01e0, 16.3098310e0, 11.68477889e0, "235.755"},
           {2, 2, 19.32e0, 8.9289802e0, 11.71529575e0, "237.555"},
           // DATA (NJ(J),MJ(J),HS(J),PHASE(J),FREQ(J),NUMARG(J),J=46,60)
           {2, 2, 1.30e0, 5.0519065e0, 11.73249771e0, "238.554"},
           {2, 2, -1.02e0, 15.8350306e0, 11.89560406e0, "244.656"},
           {2, 2, -4.51e0, 8.6624178e0, 11.91188181e0, "245.645"},
           {2, 2, 120.99e0, 11.9579569e0, 11.91280603e0, "245.655"},
           {2, 2, 1.13e0, 8.0808832e0, 11.93000800e0, "246.654"},
           {2, 2, 22.98e0, 4.5771061e0, 11.94332289e0, "247.455"},
           {2, 2, 1.06e0, 0.7000324e0, 11.96052486e0, "248.454"},
           {2, 2, -1.90e0, 14.9869335e0, 12.11031632e0, "253.755"},
           {2, 2, -2.18e0, 11.4831564e0, 12.12363121e0, "254.556"},
           {2, 2, -23.58e0, 4.3105437e0, 12.13990896e0, "255.545"},
           {2, 2, 631.92e0, 7.6060827e0, 12.14083318e0, "255.555"},
           {2, 2, 1.92e0, 3.7290090e0, 12.15803515e0, "256.554"},
           {2, 2, -4.66e0, 10.6350594e0, 12.33834347e0, "263.655"},
           {2, 2, -17.86e0, 3.2542086e0, 12.36886033e0, "265.455"},
           {2, 2, 4.47e0, 12.7336164e0, 12.37274905e0, "265.655"},
           // DATA (NJ(J),MJ(J),HS(J),PHASE(J),FREQ(J),NUMARG(J),J=61,71)
           {2, 2, 1.97e0, 16.0291555e0, 12.37367327e0, "265.665"},
           {2, 2, 17.20e0, 10.1602590e0, 12.54916865e0, "272.556"},
           {2, 2, 294.00e0, 6.2831853e0, 12.56637061e0, "273.555"},
           {2, 2, -2.46e0, 2.4061116e0, 12.58357258e0, "274.554"},
           {2, 2, -1.02e0, 5.0862033e0, 12.59985198e0, "275.545"},
           {2, 2, 79.96e0, 8.3817423e0, 12.60077620e0, "275.555"},
           {2, 2, 23.83e0, 11.6772814e0, 12.60170041e0, "275.565"},
           {2, 2, 2.59e0, 14.9728205e0, 12.60262463e0, "275.575"},
           {2, 2, 4.47e0, 4.0298682e0, 12.82880334e0, "285.455"},
           {2, 2, 1.95e0, 7.3254073e0, 12.82972756e0, "285.465"},
           {2, 2, 1.17e0, 9.1574019e0, 13.06071921e0, "295.555"}};
  int nlines = sizeof(x) / sizeof(x[0]);

  double anm[2][4][3], bnm[2][4][3];

  // Compute the time dependent potential matrix
  double dt60, pinm, alpha;
  int nn, mm;
  for (int k = -1; k < 2; k++) {
    dt60 = (dmjd - k * dt) - d1960;
    anm[0][1][k + 1] = anm[0][2][k + 1] = 0e0;
    bnm[0][1][k + 1] = bnm[0][2][k + 1] = 0e0;
    for (int j = 0; j < nlines; j++) {
      nn = x[j].nj;
      mm = x[j].mj;
      pinm = ((double)((nn + mm) % 2)) * TWOPI / 4e0;
      alpha = std::fmod(x[j].phase - pinm, TWOPI) +
              std::fmod(x[j].freq * dt60, TWOPI);
      anm[nn - 2][mm][k + 1] += x[j].hs * std::cos(alpha);
      bnm[nn - 2][mm][k + 1] -= x[j].hs * std::sin(alpha);
    }
  }

  double p[3][2], q[3][2], ap, am, bp, bm;
  // orthogonalize the response terms
  for (int m = 0; m < 2; m++) {
    ap = anm[0][m + 1][2] + anm[0][m + 1][0];
    am = anm[0][m + 1][2] - anm[0][m + 1][0];
    bp = bnm[0][m + 1][2] + bnm[0][m + 1][0];
    bm = bnm[0][m + 1][2] - bnm[0][m + 1][0];
    p[0][m] = sp[m][0] * anm[0][m + 1][1];
    p[1][m] = sp[m][1] * anm[0][m + 1][1] - sp[m][2] * ap;
    p[2][m] = sp[m][3] * anm[0][m + 1][1] - sp[m][4] * ap + sp[m][5] * bm;
    q[0][m] = sp[m][0] * bnm[0][m + 1][1];
    q[1][m] = sp[m][1] * bnm[0][m + 1][1] - sp[m][2] * bp;
    q[2][m] = sp[m][3] * bnm[0][m + 1][1] - sp[m][4] * bp - sp[m][5] * am;
    anm[0][m + 1][0] = p[0][m];
    anm[0][m + 1][1] = p[1][m];
    anm[0][m + 1][2] = p[2][m];
    bnm[0][m + 1][0] = q[0][m];
    bnm[0][m + 1][1] = q[1][m];
    bnm[0][m + 1][2] = q[2][m];
  }

  // fill partials vector
  int j = -1;
  for (int n = 2; n <= nmax; n++) {
    for (int m = 1; m <= n; m++) {
      for (int k = -1; k <= 1; k++) {
        h[++j] = anm[n - 2][m][k + 1];
        h[++j] = bnm[n - 2][m][k + 1];
      }
    }
  }

  // A final check !
  assert(j == 11);

  // Finished
  return 0;
}
