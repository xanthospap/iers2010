#include "fundarg.hpp"

double *dso::iers2010::fundarg(double t, double *fargs) noexcept {
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
  return fargs;
}

double *dso::iers2010::fundarg_derivs(double t, double *fargs) noexcept {
  /*  Compute the fundamental argument L */
  fargs[0] = iers2010::dfal03(t);

  /* Compute the fundamental argument LP */
  fargs[1] = iers2010::dfalp03(t);

  /* Compute the fundamental argument F. */
  fargs[2] = iers2010::dfaf03(t);

  /* Compute the fundamental argument D. */
  fargs[3] = iers2010::dfad03(t);

  /* Compute the fundamental argument OM. */
  fargs[4] = iers2010::dfaom03(t);
  
  /* Finished. */
  return fargs;
}
