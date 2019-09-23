#include "iers2010.hpp"
#ifdef USE_EXTERNAL_CONSTS
#include "gencon.hpp"
#endif

/// @details This subroutine computes the lunisolar fundamental arguments.
///          The model used is from Simon et al. (1994) as recommended by the 
///          IERS Conventions (2010).  Refer to IERS Conventions (2010) 
///          Chapter 5 Sections 5.7.1 - 5.7.2 (pp. 57 - 59).
///          This function is a translation/wrapper for the fortran FUNDARG
///          subroutine, found here : 
///          http://maia.usno.navy.mil/conv2010/software.html
/// 
/// @param[in]  t     TT, Julian centuries since J2000 (Note 1)
/// @param[out] fargs A 5-element array containing the computed fundamental
///                   arguments, in the following order:
///                   fargs[0] -> l  : Mean anomaly of the Moon (Note 2)
///                   fargs[1] -> lp : Mean anomaly of the Sun (Note 2)
///                   fargs[2] -> f  : L - OM (Notes 2 and 3)
///                   fargs[3] -> d  : Mean elongation of the Moon from the Sun
///                                    (Note 2)
///                   fargs[4] -> om : Mean longitude of the ascending node of
///                                    the Moon (Note 2)
/// @return           An integer value always 0.
/// 
/// @note
///       -# Though T is strictly TDB, it is usually more convenient to use
///          TT, which makes no significant difference.  Julian centuries since
///          J2000 is (JD - 2451545.0)/36525.
///       -# The expression used is as adopted in IERS Conventions (2010) and
///          is from Simon et al. (1994).  Arguments are in radians.
///       -# L in this instance is the Mean Longitude of the Moon. OM is the 
///          Mean longitude of the ascending node of the Moon.
///       -# Status: Canonical model
/// 
/// @version 25.02.2010
/// 
/// @cite iers2010
/// @cite simon94
int
iers2010::fundarg(double t, double* fargs)
{

// Set constants
#ifndef USE_EXTERNAL_CONSTS
  // Arcseconds to radians
  constexpr double DAS2R   (4.848136811095359935899141e-6);
  // Arcseconds in a full circle
  constexpr double TURNAS  (1296000e0);
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  //  Compute the fundamental argument L.
  fargs[0] = std::fmod(485868.249036e0 +
          t*( 1717915923.2178e0 +
          t*(         31.8792e0 +
          t*(          0.051635e0 +
          t*(        - 0.00024470e0 )))), TURNAS ) * DAS2R;

  // Compute the fundamental argument LP.
  fargs[1] = std::fmod(1287104.79305e0 +
              t*( 129596581.0481e0 +
              t*(       - 0.5532e0 +
              t*(         0.000136e0 +
              t*(       - 0.00001149e0 )))), TURNAS ) * DAS2R;

  // Compute the fundamental argument F.
  fargs[2] = std::fmod(335779.526232e0 +
              t*( 1739527262.8478e0 +
              t*(       - 12.7512e0 +
              t*(       -  0.001037e0 +
              t*(          0.00000417e0 )))), TURNAS ) * DAS2R;

  // Compute the fundamental argument D.
  fargs[3] = std::fmod(1072260.70369e0 +
              t*( 1602961601.2090e0 +
              t*(        - 6.3706e0 +
              t*(          0.006593e0 +
              t*(        - 0.00003169e0 )))), TURNAS ) * DAS2R;

  // Compute the fundamental argument OM.
  fargs[4] = std::fmod(450160.398036e0 +
              t*( - 6962890.5431e0 +
              t*(         7.4722e0 +
              t*(         0.007702e0 +
              t*(       - 0.00005939e0 )))), TURNAS ) * DAS2R;

  // Finished.
  return 0;
}
