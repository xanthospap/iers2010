#include "iersc.hpp"
#include "iers2010.hpp"
#include <cmath>

//constexpr const double halfpi = iers2010::DPI / 2e0;
constexpr const double secrad = iers2010::DPI / (180e0*3600e0);
constexpr const double secfac = 1296000e0; // see 5.4 in IERS Conventions 2010
constexpr const double dayscentury = 36525e0;

/// @brief Tables 8.2a and 8.2b from IERS-2010 (Doodson numbers and periods
///        not relevant here).
/// @see IERS Conventions 8.2 and the routine interp.f at 
///      https://hpiers.obspm.fr/iers/models/interp.readme &
///      https://hpiers.obspm.fr/iers/models/interp.f (SUBROUTINE PMUT1_OCEANS)
struct {
  ///< Multipliers of GMST+pi and Delaunay arguments.
  int narg[6];
  ///< Oceanic tidal terms present in x (microas),y(microas),ut1(microseconds)
  double xsin,xcos,ysin,ycos,utsin,utcos;
} data[] = {
   {{1,-1, 0,-2,-2,-2},  -0.05e0,   0.94e0,  -0.94e0,  -0.05e0,  0.396e0, -0.078e0},
   {{1,-2, 0,-2, 0,-1},   0.06e0,   0.64e0,  -0.64e0,   0.06e0,  0.195e0, -0.059e0},
   {{1,-2, 0,-2, 0,-2},   0.30e0,   3.42e0,  -3.42e0,   0.30e0,  1.034e0, -0.314e0},
   {{1, 0, 0,-2,-2,-1},   0.08e0,   0.78e0,  -0.78e0,   0.08e0,  0.224e0, -0.073e0},
   {{1, 0, 0,-2,-2,-2},   0.46e0,   4.15e0,  -4.15e0,   0.45e0,  1.187e0, -0.387e0},
   {{1,-1, 0,-2, 0,-1},   1.19e0,   4.96e0,  -4.96e0,   1.19e0,  0.966e0, -0.474e0},
   {{1,-1, 0,-2, 0,-2},   6.24e0,  26.31e0, -26.31e0,   6.23e0,  5.118e0, -2.499e0},
   {{1, 1, 0,-2,-2,-1},   0.24e0,   0.94e0,  -0.94e0,   0.24e0,  0.172e0, -0.090e0},
   {{1, 1, 0,-2,-2,-2},   1.28e0,   4.99e0,  -4.99e0,   1.28e0,  0.911e0, -0.475e0},
   {{1, 0, 0,-2, 0, 0},  -0.28e0,  -0.77e0,   0.77e0,  -0.28e0, -0.093e0,  0.070e0},
   {{1, 0, 0,-2, 0,-1},   9.22e0,  25.06e0, -25.06e0,   9.22e0,  3.025e0, -2.280e0},
   {{1, 0, 0,-2, 0,-2},  48.82e0, 132.91e0,-132.90e0,  48.82e0, 16.020e0,-12.069e0},
   {{1,-2, 0, 0, 0, 0},  -0.32e0,  -0.86e0,   0.86e0,  -0.32e0, -0.103e0,  0.078e0},
   {{1, 0, 0, 0,-2, 0},  -0.66e0,  -1.72e0,   1.72e0,  -0.66e0, -0.194e0,  0.154e0},
   {{1,-1, 0,-2, 2,-2},  -0.42e0,  -0.92e0,   0.92e0,  -0.42e0, -0.083e0,  0.074e0},
   {{1, 1, 0,-2, 0,-1},  -0.30e0,  -0.64e0,   0.64e0,  -0.30e0, -0.057e0,  0.050e0},
   {{1, 1, 0,-2, 0,-2},  -1.61e0,  -3.46e0,   3.46e0,  -1.61e0, -0.308e0,  0.271e0},
   {{1,-1, 0, 0, 0, 0},  -4.48e0,  -9.61e0,   9.61e0,  -4.48e0, -0.856e0,  0.751e0},
   {{1,-1, 0, 0, 0,-1},  -0.90e0,  -1.93e0,   1.93e0,  -0.90e0, -0.172e0,  0.151e0},
   {{1, 1, 0, 0,-2, 0},  -0.86e0,  -1.81e0,   1.81e0,  -0.86e0, -0.161e0,  0.137e0},
   {{1, 0,-1,-2, 2,-2},   1.54e0,   3.03e0,  -3.03e0,   1.54e0,  0.315e0, -0.189e0},
   {{1, 0, 0,-2, 2,-1},  -0.29e0,  -0.58e0,   0.58e0,  -0.29e0, -0.062e0,  0.035e0},
   {{1, 0, 0,-2, 2,-2},  26.13e0,  51.25e0, -51.25e0,  26.13e0,  5.512e0, -3.095e0},
   {{1, 0, 1,-2, 2,-2},  -0.22e0,  -0.42e0,   0.42e0,  -0.22e0, -0.047e0,  0.025e0},
   {{1, 0,-1, 0, 0, 0},  -0.61e0,  -1.20e0,   1.20e0,  -0.61e0, -0.134e0,  0.070e0},
   {{1, 0, 0, 0, 0, 1},   1.54e0,   3.00e0,  -3.00e0,   1.54e0,  0.348e0, -0.171e0},
   {{1, 0, 0, 0, 0, 0}, -77.48e0,-151.74e0, 151.74e0, -77.48e0,-17.620e0,  8.548e0},
   {{1, 0, 0, 0, 0,-1}, -10.52e0, -20.56e0,  20.56e0, -10.52e0, -2.392e0,  1.159e0},
   {{1, 0, 0, 0, 0,-2},   0.23e0,   0.44e0,  -0.44e0,   0.23e0,  0.052e0, -0.025e0},
   {{1, 0, 1, 0, 0, 0},  -0.61e0,  -1.19e0,   1.19e0,  -0.61e0, -0.144e0,  0.065e0},
   {{1, 0, 0, 2,-2, 2},  -1.09e0,  -2.11e0,   2.11e0,  -1.09e0, -0.267e0,  0.111e0},
   {{1,-1, 0, 0, 2, 0},  -0.69e0,  -1.43e0,   1.43e0,  -0.69e0, -0.288e0,  0.043e0},
   {{1, 1, 0, 0, 0, 0},  -3.46e0,  -7.28e0,   7.28e0,  -3.46e0, -1.610e0,  0.187e0},
   {{1, 1, 0, 0, 0,-1},  -0.69e0,  -1.44e0,   1.44e0,  -0.69e0, -0.320e0,  0.037e0},
   {{1, 0, 0, 0, 2, 0},  -0.37e0,  -1.06e0,   1.06e0,  -0.37e0, -0.407e0, -0.005e0},
   {{1, 2, 0, 0, 0, 0},  -0.17e0,  -0.51e0,   0.51e0,  -0.17e0, -0.213e0, -0.005e0},
   {{1, 0, 0, 2, 0, 2},  -1.10e0,  -3.42e0,   3.42e0,  -1.09e0, -1.436e0, -0.037e0},
   {{1, 0, 0, 2, 0, 1},  -0.70e0,  -2.19e0,   2.19e0,  -0.70e0, -0.921e0, -0.023e0},
   {{1, 0, 0, 2, 0, 0},  -0.15e0,  -0.46e0,   0.46e0,  -0.15e0, -0.193e0, -0.005e0},
   {{1, 1, 0, 2, 0, 2},  -0.03e0,  -0.59e0,   0.59e0,  -0.03e0, -0.396e0, -0.024e0},
   {{1, 1, 0, 2, 0, 1},  -0.02e0,  -0.38e0,   0.38e0,  -0.02e0, -0.253e0, -0.015e0},
   {{2,-3, 0,-2, 0,-2},  -0.49e0,  -0.04e0,   0.63e0,   0.24e0, -0.089e0, -0.011e0},
   {{2,-1, 0,-2,-2,-2},  -1.33e0,  -0.17e0,   1.53e0,   0.68e0, -0.224e0, -0.032e0},
   {{2,-2, 0,-2, 0,-2},  -6.08e0,  -1.61e0,   3.13e0,   3.35e0, -0.637e0, -0.177e0},
   {{2, 0, 0,-2,-2,-2},  -7.59e0,  -2.05e0,   3.44e0,   4.23e0, -0.745e0, -0.222e0},
   {{2, 0, 1,-2,-2,-2},  -0.52e0,  -0.14e0,   0.22e0,   0.29e0, -0.049e0, -0.015e0},
   {{2,-1,-1,-2, 0,-2},   0.47e0,   0.11e0,  -0.10e0,  -0.27e0,  0.033e0,  0.013e0},
   {{2,-1, 0,-2, 0,-1},   2.12e0,   0.49e0,  -0.41e0,  -1.23e0,  0.141e0,  0.058e0},
   {{2,-1, 0,-2, 0,-2}, -56.87e0, -12.93e0,  11.15e0,  32.88e0, -3.795e0, -1.556e0},
   {{2,-1, 1,-2, 0,-2},  -0.54e0,  -0.12e0,   0.10e0,   0.31e0, -0.035e0, -0.015e0},
   {{2, 1, 0,-2,-2,-2}, -11.01e0,  -2.40e0,   1.89e0,   6.41e0, -0.698e0, -0.298e0},
   {{2, 1, 1,-2,-2,-2},  -0.51e0,  -0.11e0,   0.08e0,   0.30e0, -0.032e0, -0.014e0},
   {{2,-2, 0,-2, 2,-2},   0.98e0,   0.11e0,  -0.11e0,  -0.58e0,  0.050e0,  0.022e0},
   {{2, 0,-1,-2, 0,-2},   1.13e0,   0.11e0,  -0.13e0,  -0.67e0,  0.056e0,  0.025e0},
   {{2, 0, 0,-2, 0,-1},  12.32e0,   1.00e0,  -1.41e0,  -7.31e0,  0.605e0,  0.266e0},
   {{2, 0, 0,-2, 0,-2},-330.15e0, -26.96e0,  37.58e0, 195.92e0,-16.195e0, -7.140e0},
   {{2, 0, 1,-2, 0,-2},  -1.01e0,  -0.07e0,   0.11e0,   0.60e0, -0.049e0, -0.021e0},
   {{2,-1, 0,-2, 2,-2},   2.47e0,  -0.28e0,  -0.44e0,  -1.48e0,  0.111e0,  0.034e0},
   {{2, 1, 0,-2, 0,-2},   9.40e0,  -1.44e0,  -1.88e0,  -5.65e0,  0.425e0,  0.117e0},
   {{2,-1, 0, 0, 0, 0},  -2.35e0,   0.37e0,   0.47e0,   1.41e0, -0.106e0, -0.029e0},
   {{2,-1, 0, 0, 0,-1},  -1.04e0,   0.17e0,   0.21e0,   0.62e0, -0.047e0, -0.013e0},
   {{2, 0,-1,-2, 2,-2},  -8.51e0,   3.50e0,   3.29e0,   5.11e0, -0.437e0, -0.019e0},
   {{2, 0, 0,-2, 2,-2},-144.13e0,  63.56e0,  59.23e0,  86.56e0, -7.547e0, -0.159e0},
   {{2, 0, 1,-2, 2,-2},   1.19e0,  -0.56e0,  -0.52e0,  -0.72e0,  0.064e0,  0.000e0},
   {{2, 0, 0, 0, 0, 1},   0.49e0,  -0.25e0,  -0.23e0,  -0.29e0,  0.027e0, -0.001e0},
   {{2, 0, 0, 0, 0, 0}, -38.48e0,  19.14e0,  17.72e0,  23.11e0, -2.104e0,  0.041e0},
   {{2, 0, 0, 0, 0,-1}, -11.44e0,   5.75e0,   5.32e0,   6.87e0, -0.627e0,  0.015e0},
   {{2, 0, 0, 0, 0,-2},  -1.24e0,   0.63e0,   0.58e0,   0.75e0, -0.068e0,  0.002e0},
   {{2, 1, 0, 0, 0, 0},  -1.77e0,   1.79e0,   1.71e0,   1.04e0, -0.146e0,  0.037e0},
   {{2, 1, 0, 0, 0,-1},  -0.77e0,   0.78e0,   0.75e0,   0.45e0, -0.064e0,  0.017e0},
   {{2, 0, 0, 2, 0, 2},  -0.33e0,   0.62e0,   0.65e0,   0.19e0, -0.049e0,  0.018e0}
};

constexpr const int nl = sizeof(data) / sizeof(data[0]);
static_assert(nl == 71);

int iers2010::interp::pmut1_oceans(double t /*in julian centuries*/,
                                   double &cor_x, double &cor_y,
                                   double &cor_ut1, double &cor_lod) noexcept {
  double arg[6];  // Array of the tidal arguments
  double darg[6]; // Array of their time derivative
  const double t2 = t * t;
  const double t3 = t2 * t;
  const double t4 = t3 * t;

  // Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
  arg[0] = (67310.54841e0 + (876600e0 * 3600e0 + 8640184.812866e0) * t +
            0.093104e0 * t2 - 6.2e0 - 6e0 * t3) *
               15.0e0 +
           648000.0e0;
  arg[0] = std::fmod(arg[0], secfac) * secrad;

  arg[1] = -0.00024470e0 * t4 + 0.051635e0 * t3 + 31.8792e0 * t2 +
           1717915923.2178e0 * t + 485868.249036e0;
  arg[1] = std::fmod(arg[1], secfac) * secrad;

  arg[2] = -0.00001149e0 * t4 - 0.000136e0 * t3 - 0.5532e0 * t2 +
           129596581.0481e0 * t + 1287104.79305e0;
  arg[2] = std::fmod(arg[2], secfac) * secrad;

  arg[3] = 0.00000417e0 * t4 - 0.001037e0 * t3 - 12.7512e0 * t2 +
           1739527262.8478e0 * t + 335779.526232e0;
  arg[3] = std::fmod(arg[3], secfac) * secrad;

  arg[4] = -0.00003169e0 * t4 + 0.006593e0 * t3 - 6.3706e0 * t2 +
           1602961601.2090e0 * t + 1072260.70369e0;
  arg[4] = std::fmod(arg[4], secfac) * secrad;

  arg[5] = -0.00005939e0 * t4 + 0.007702e0 * t3 + 7.4722e0 * t2 -
           6962890.2665e0 * t + 450160.398036e0;
  arg[5] = std::fmod(arg[5], secfac) * secrad;

  darg[0] = (876600e0 * 3600e0 + 8640184.812866e0 + 2.e0 * 0.093104e0 * t -
             3.e0 * 6.2e0 - 6 * t2) *
            15.e0;
  darg[0] = darg[0] * secrad / dayscentury; // rad/day

  darg[1] = -4.e0 * 0.00024470e0 * t3 + 3.e0 * 0.051635e0 * t2 +
            2.e0 * 31.8792e0 * t + 1717915923.2178e0;
  darg[1] = darg[1] * secrad / dayscentury; // rad/day

  darg[2] = -4.e0 * 0.00001149e0 * t3 - 3.e0 * 0.000136e0 * t2 -
            2.e0 * 0.5532e0 * t + 129596581.0481e0;
  darg[2] = darg[2] * secrad / dayscentury; // rad/day

  darg[3] = 4.e0 * 0.00000417e0 * t3 - 3.e0 * 0.001037e0 * t2 -
            2.e0 * 12.7512e0 * t + 1739527262.8478e0;
  darg[3] = darg[3] * secrad / dayscentury; // rad/day

  darg[4] = -4.e0 * 0.00003169e0 * t3 + 3.e0 * 0.006593e0 * t2 -
            2.e0 * 6.3706e0 * t + 1602961601.2090e0;
  darg[4] = darg[4] * secrad / dayscentury; // rad/day

  darg[5] = -4.e0 * 0.00005939e0 * t3 + 3.e0 * 0.007702e0 * t2 +
            2.e0 * 7.4722e0 * t - 6962890.2665e0;
  darg[5] = darg[5] * secrad / dayscentury; // rad/day

  // Corrections
  cor_x = 0e0;
  cor_y = 0e0;
  cor_ut1 = 0e0;
  cor_lod = 0e0;

  for (int j = 0; j < nl; j++) {
    double ag = 0e0;
    double dag = 0e0;
    for (int i = 0; i < 6; i++) {
      ag += static_cast<double>(data[j].narg[i]) * arg[i];
      dag += static_cast<double>(data[j].narg[i]) * darg[i];
    }
    ag = std::fmod(ag, iers2010::D2PI);
    const double cag = std::cos(ag);
    const double sag = std::sin(ag);

    cor_x += data[j].xcos * cag + data[j].xsin * sag;
    cor_y += data[j].ycos * cag + data[j].ysin * sag;
    cor_ut1 += data[j].utcos * cag + data[j].utsin * sag;
    cor_lod -= (-data[j].utcos * sag + data[j].utsin * cag) * dag;
  }

  // Units
  cor_x *= 1e-6;   // arcseconds ('')
  cor_y *= 1e-6;   // arcseconds ('')
  cor_ut1 *= 1e-6; // seconds (s)
  cor_lod *= 1e-6; // seconds (s)

  return 0;
}
