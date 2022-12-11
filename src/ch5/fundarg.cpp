#include "iers2010.hpp"

namespace {
constexpr const double fTURNAS = 1296000e0;
constexpr const double fDAS2R = 4.848136811095359935899141e-6;
}

int iers2010::fundarg(double t, double *fargs) noexcept {

  //  Compute the fundamental argument L.
  fargs[0] = std::fmod(485868.249036e0 +
                           t * (1717915923.2178e0 +
                                t * (31.8792e0 +
                                     t * (0.051635e0 + t * (-0.00024470e0)))),
                       fTURNAS) *
             fDAS2R;

  // Compute the fundamental argument LP.
  fargs[1] = std::fmod(1287104.79305e0 +
                           t * (129596581.0481e0 +
                                t * (-0.5532e0 +
                                     t * (0.000136e0 + t * (-0.00001149e0)))),
                       fTURNAS) *
             fDAS2R;

  // Compute the fundamental argument F.
  fargs[2] = std::fmod(335779.526232e0 +
                           t * (1739527262.8478e0 +
                                t * (-12.7512e0 +
                                     t * (-0.001037e0 + t * (0.00000417e0)))),
                       fTURNAS) *
             fDAS2R;

  // Compute the fundamental argument D.
  fargs[3] = std::fmod(1072260.70369e0 +
                           t * (1602961601.2090e0 +
                                t * (-6.3706e0 +
                                     t * (0.006593e0 + t * (-0.00003169e0)))),
                       fTURNAS) *
             fDAS2R;

  // Compute the fundamental argument OM.
  fargs[4] = std::fmod(450160.398036e0 +
                           t * (-6962890.5431e0 +
                                t * (7.4722e0 +
                                     t * (0.007702e0 + t * (-0.00005939e0)))),
                       fTURNAS) *
             fDAS2R;

  // Finished.
  return 0;
}
