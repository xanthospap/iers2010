#include "fundarg.hpp"

void iers2010::details::fundarg(double t, double *fargs) noexcept {
  /*  Compute the fundamental argument L */
  fargs[0] = iers2010::fal03(t);

  /* Compute the fundamental argument LP */
  fargs[1] = iers2010::falp03(t);

  /* Compute the fundamental argument F. */
  fargs[2] = iers2010::faf03(t);

  /* Compute the fundamental argument D. */
  fargs[3] = iers2010::fad03(t);

  /* Compute the fundamental argument OM. */
  fargs[4] = iers2010::faom03(t);

  /* Finished. */
  return;
}

void iers2010::details::fundarg_derivs(double t, double *fargs) noexcept {
  // Note that obviously we are not normalizing angles here! (aka no fmod,
  // anp, etc ...)

  //  Compute the derivative of fundamental argument L.
  fargs[0] = (1717915923.2178e0 +
              t * (31.8792e0 * 2e0 +
                   t * (0.051635e0 * 3e0 + t * (-0.00024470e0 * 4e0)))) *
             iers2010::DAS2R;

  // Compute the derivative of fundamental argument LP.
  fargs[1] = (129596581.0481e0 +
              t * (-0.5532e0 * 2e0 +
                   t * (0.000136e0 * 3e0 + t * (-0.00001149e0 * 4e0)))) *
             iers2010::DAS2R;

  // Compute the derivative of fundamental argument F.
  fargs[2] = (1739527262.8478e0 +
              t * (-12.7512e0 * 2e0 +
                   t * (-0.001037e0 * 3e0 + t * (0.00000417e0 * 4e0)))) *
             iers2010::DAS2R;

  // Compute the derivative of fundamental argument D.
  fargs[3] = (1602961601.2090e0 +
              t * (-6.3706e0 * 2e0 +
                   t * (0.006593e0 * 3e0 + t * (-0.00003169e0 * 4e0)))) *
             iers2010::DAS2R;

  // Compute the derivative of fundamental argument OM.
  fargs[4] = (-6962890.5431e0 +
              t * (7.4722e0 * 2e0 +
                   t * (0.007702e0 * 3e0 + t * (-0.00005939e0 * 4e0)))) *
             iers2010::DAS2R;

  // Finished.
  return;
}
