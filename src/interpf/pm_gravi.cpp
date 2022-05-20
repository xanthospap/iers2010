#include "iersc.hpp"
#include "iers2010.hpp"
#include <cmath>

// constexpr const double halfpi = iers2010::DPI / 2e0;
constexpr const double secrad = iers2010::DPI / (180e0*3600e0);
constexpr const double secfac = 1296000e0; // see 5.4 in IERS Conventions 2010

struct {
  // Multipliers of GMST+pi and Delaunay arguments.
  int narg[6];
  // Diurnal lunisolar tidal terms present in x (microas),y(microas)
  double xsin,xcos,ysin,ycos;
} data[] = {
     {{ 1,-1, 0,-2, 0,-1},    -.44e0,   .25e0,   -.25e0,  -.44e0},
     {{ 1,-1, 0,-2, 0,-2},   -2.31e0,  1.32e0,  -1.32e0, -2.31e0},
     {{ 1, 1, 0,-2,-2,-2},    -.44e0,   .25e0,   -.25e0,  -.44e0},
     {{ 1, 0, 0,-2, 0,-1},   -2.14e0,  1.23e0,  -1.23e0, -2.14e0},
     {{ 1, 0, 0,-2, 0,-2},  -11.36e0,  6.52e0,  -6.52e0,-11.36e0},
     {{ 1,-1, 0, 0, 0, 0},     .84e0,  -.48e0,    .48e0,   .84e0},
     {{ 1, 0, 0,-2, 2,-2},   -4.76e0,  2.73e0,  -2.73e0, -4.76e0},
     {{ 1, 0, 0, 0, 0, 0},   14.27e0, -8.19e0,   8.19e0, 14.27e0},
     {{ 1, 0, 0, 0, 0,-1},    1.93e0, -1.11e0,   1.11e0,  1.93e0},
     {{ 1, 1, 0, 0, 0, 0},     .76e0,  -.43e0,    .43e0,   .76e0}
};

constexpr const int nl = sizeof(data) / sizeof(data[0]);
static_assert(nl == 10);

int iers2010::interp::pm_gravi(double t, double &cor_x,
                               double &cor_y) noexcept {
  double arg[6]; // Array of the tidal arguments
  const double t2 = t * t;
  const double t3 = t2 * t;
  const double t4 = t3 * t;

  // Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
  arg[0] = (67310.54841e0 + (876600e0 * 3600e0 + 8640184.812866e0) * t +
            0.093104e0 * t2 - 6.2e0 - 6 * t3) *
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

  // corrections
  cor_x = 0e0;
  cor_y = 0e0;
  for (int j = 0; j < nl; j++) {
    double ag = 0e0;
    for (int i = 0; i < 6; i++)
      ag += static_cast<double>(data[j].narg[i]) * arg[i];
    ag = std::fmod(ag, iers2010::D2PI);

    const double cag = std::cos(ag);
    const double sag = std::sin(ag);
    cor_x += data[j].xcos * cag + data[j].xsin * sag;
    cor_y += data[j].ycos * cag + data[j].ysin * sag;
  }

  // Units
  // cor_x *= 1e-6; // arcseconds ('')
  // cor_y *= 1e-6; // arcseconds ('')
  return 0;
}
