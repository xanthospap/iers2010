#include "iers2010.hpp"

#ifdef USE_EXTERNAL_CONSTS
  #include "gencon.hpp"
#endif

#ifdef QUICK_EXIT
  #include <algorithm>
#endif

/**
 * @details The purpose of the function is to compute the angular astronomical
 *          argument, which depends on time, for 11 tidal argument calculations. 
 *          The order of the 11 angular quantities in vector angle are given
 *          below:
 *          01-M2, 02-S2, 03-N2, 04-K2, 05-K1, 06-O1, 07-P1, 08-Q1, 09-Mf,
 *          10-Mm, 11-Ssa (See Reference 1) 
 *          This function is a translation/wrapper for the fortran ARG2
 *          subroutine, found here : 
 *          http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  iyear Four digit year (Note 1)
 * @param[in]  day   Day of Year Greenwich Time (Note 2)
 * @param[out] angle Angular argument for Schwiderski computation, in radians 
 *                   (Notes 3, 4 and 5)
 * @return           An integer value which can be:
 *                   Returned Value | Status
 *                   ---------------|-------------------------------------------
 *                               -1 | Error; Invalid year
 *                                0 | All ok; a new value has been computed
 * (Only when QUICK_EXIT enabled) 1 | All ok; previous value for angle used.
 *
 * @note 
 *  -# This subroutine is valid only after 1973 CE.  A validation
 *     test has been added to stop the subroutine if an invalid
 *     year is used.
 *  -# Example: 32.5 for February 1 12 Noon
 *     Example: 1.25 for January 1 6 AM
 *  -# Ocean loading phases computed from Schwiderski's models
 *     refer to the phase of the associated solid Earth tide 
 *     generating potential at the zero meridian according to <br>
 *      OL_DR = OL_AMP ' COS (SE_PHASE" - OL_PHASE) <br>
 *     where OL = OCEAN LOADING TIDE,<br>
 *           SE = SOLID EARTH TIDE GENERATING POTENTIAL.<br>
 *     If the harmonic tide development of Cartwright, et al. 
 *     (CTE) (1971, 1973) is used, make sure that SE_PHASE"
 *     take into account: 
 *     - the sign of SE_AMP in the tables of Cartwright et al.
 *     - that CTE'S SE_PHASE refers to a sine rather than a 
 *       cosine function if (N+M) = (DEGREE + ORDER) of the tide
 *       spherical harmonic is odd.
 *     i.e. SE_PHASE" = TAU(T) ' N1 + S(T) ' N2 + H(T) ' N3 <br>
 *                 + P(T) ' N4 + N'(T) ' N5 + PS(T) ' N6 <br>
 *                 + PI   If CTE'S amplitude coefficient < 0 <br>
 *                 + PI/2 If (DEGREE + N1) is odd <br>
 *     where TAU ... PS = astronomical arguments,<br>
 *           N1 ... N6 = CTE'S argument numbers.<br>
 *     Most tide generating software compute SE_PHASE" (for use
 *     with cosines).
 *  -# The double precision change from the original routine ARG.f
 *     to ARG2.F yields output differences on the order of 10^-9 radians.
 *  -# The input array \c angle must be able to hold at least 11 doubles.
 *  -# Status: Canonical model
 *
 *  Test case:
 *     given input: IYEAR = 2008 
 *                  DAY = 311.5 (November 6 Noon)
 *     expected output: ANGLE(1)  = 2.849663065753787805D0  rad
 *                      ANGLE(2)  = 6.28318080000000023D0   rad
 *                      ANGLE(3)  = 4.926040134021299366D0  rad
 *                      ANGLE(4)  = 1.608450491115348768D0  rad
 *                      ANGLE(5)  = 2.375021572352622456D0  rad
 *                      ANGLE(6)  = 0.4746414933980958040D0 rad
 *                      ANGLE(7)  = 3.908159227647345801D0  rad
 *                      ANGLE(8)  = 2.551018561669245344D0  rad
 *                      ANGLE(9)  = 5.041990012540757959D0  rad 
 *                      ANGLE(10) = 4.206816878908014701D0  rad 
 *                      ANGLE(11) = 1.608463638294885811D0  rad 
 *
 * @version  2011 October   07
 *
 * @cite iers2010
 * Schwiderski, E., 1983, "Atlas of Ocean Tidal Charts and Maps, Part I:
 * The Semidiurnal Principal Lunar Tide M2," Marine Geodesy, 6, pp. 219-256.
 *
 */
int iers2010::arg2 (const int& iyear,const double& day,double* angle) 
{
  // check for quick exit
  #ifdef QUICK_EXIT
    static int    iyear_previous = .0e0;
    static double day_previous = -999.0;
    static double angle_previous[11];
    if (iyear==iyear_previous) {
      if (fabs(day_previous-day)<DATE_MAX_DIFF) {
        std::copy (angle_previous,angle_previous+11,angle);
        return 1;
      }
    }
  #endif

  // Constants
  const int k = 11;
  const int iymin = 1974;
  
  #ifdef USE_EXTERNAL_CONSTS
    constexpr double TWOPI   (D2PI);
  #else
    constexpr double TWOPI   (6.283185307179586476925287e0);
  #endif

  /*  ----------------------------------------------
   *  Speed of all terms given in radians per second
   *  ---------------------------------------------- */
  static const double speed[] = {
    1.40519e-4,
    1.45444e-4,
    1.37880e-4,
    1.45842e-4,
    0.72921e-4,
    0.67598e-4,
    0.72523e-4,
    0.64959e-4,
    0.053234e-4,
    0.026392e-4,
    0.003982e-4,
  };

  /* these are not used for now ...
  const double sigm2  = 1.40519e-4;
  const double sigs2  = 1.45444e-4;
  const double sign2  = 1.37880e-4;
  const double sigk2  = 1.45842e-4;
  const double sigk1  = 0.72921e-4;
  const double sigo1  = 0.67598e-4;
  const double sigp1  = 0.72523e-4;
  const double sigq1  = 0.64959e-4;
  const double sigmf  = 0.053234e-4;
  const double sigmm  = 0.026392e-4;
  const double sigssa = 0.003982e-4;
  */
 
  static const double angfac[][4] = {
    { 0.200e+01, -0.200e+01,  0.000e+00,  0.000e+00 }, 
    { 0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00 }, 
    { 0.200e+01, -0.300e+01,  0.100e+01,  0.000e+00 }, 
    { 0.200e+01,  0.000e+00,  0.000e+00,  0.000e+00 }, 
    { 0.100e+01,  0.000e+00,  0.000e+00,  0.250e+00 }, 
    { 0.100e+01, -0.200e+01,  0.000e+00, -0.250e+00 }, 
    {-0.100e+01,  0.000e+00,  0.000e+00, -0.250e+00 }, 
    { 0.100e+01, -0.300e+01,  0.100e+01, -0.250e+00 }, 
    { 0.000e+00,  0.200e+01,  0.000e+00,  0.000e+00 }, 
    { 0.000e+00,  0.100e+01, -0.100e+01,  0.000e+00 }, 
    { 0.200e+01,  0.000e+00,  0.000e+00,  0.000e+00 }
  };

  const double dtr = 0.174532925199e-1;

  //  Validate year
  if (iyear<iymin) return -1;

  // Initialize day of year
  double id,fraction;
  fraction = modf(day,&id);

  /* ------------------------------------------
   *  Compute fractional part of day in seconds 
   * ------------------------------------------ */
  double fday = fraction * 86400e0;
  // Revision 07 October 2011: ICAPD modified 
  int icapd = (int)id + 365 * (iyear-1975) + ( (iyear-1973) /4 );
  double capt = (27392.500528e0 + 1.000000035e0 * (double)icapd) / 36525e0;
  
  /* --------------------------------------------------
   *  Compute mean longitude of Sun at beginning of day
   * -------------------------------------------------- */
  double h0 = ( 279.69668e0 + (36000.768930485e0 + 3.03e-4 * capt ) * capt ) 
    * dtr;
  
  /* ---------------------------------------------------
   *  Compute mean longitude of Moon at beginning of day 
   * --------------------------------------------------- */
  double s0 = ( ( (1.9e-6*capt-.001133e0) * capt + 481267.88314137e0) 
      * capt +270.434358e0 ) * dtr;
  
  /* ------------------------------------------------------------
   *  Compute mean longitude of lunar perigee at beginning of day 
   * ------------------------------------------------------------ */
  double p0 = ( ( (-1.2e-5*capt-.010325e0) * capt + 4069.0340329577e0 ) 
      * capt + 334.329653e0 ) * dtr;

  // Compute the tidal angle arguments
  for (int i=0;i<k;i++) {
    angle[i] = speed[i]*fday +
               angfac[i][0]*h0 +
               angfac[i][1]*s0 +
               angfac[i][2]*p0 +
               angfac[i][3]*TWOPI;
    angle[i] = fmod (angle[i],TWOPI);
    while (angle[i]<0e0) angle[i] += TWOPI;
  }
  
  // update quick exit
  #ifdef QUICK_EXIT
    iyear_previous = iyear;
    day_previous = day;
    std::copy(angle,angle+11,angle_previous);
  #endif

  // Finished.
  return 0;
}
