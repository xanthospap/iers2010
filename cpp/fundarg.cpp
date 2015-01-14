#include "iers2010.hpp"
#ifdef USE_EXTERNAL_CONSTS
  #include "gencon.hpp"
#endif

inline double _MOD_ (const double& a,const double& p)
{
  return a - ( (int)(a/p) * p );
}

/**
 * @details This subroutine computes the lunisolar fundamental arguments.
 *          The model used is from Simon et al. (1994) as recommended by the IERS
 *          Conventions (2010).  Refer to IERS Conventions (2010) Chapter 5 
 *          Sections 5.7.1 - 5.7.2 (pp. 57 - 59).
 *          This function is a translation/wrapper for the fortran FUNDARG
 *          subroutine, found here : http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  t    TT, Julian centuries since J2000 (Note 1)
 * @param[out] l    Mean anomaly of the Moon (Note 2)
 * @param[out] lp   Mean anomaly of the Sun (Note 2)
 * @param[out] f    L - OM (Notes 2 and 3)
 * @param[out] d    Mean elongation of the Moon from the Sun (Note 2)
 * @param[out] om   Mean longitude of the ascending node of the 
 *                  Moon (Note 2)
 * @return          An integer value which can be either 0 (to denote
 *                  that a new value for pm has been computed) or 1 
 *                  (to denote quick return). The returned integer has
 *                  nothing to do with the status of the function and
 *                  has no meaning if the <b>QUICK_EXIT</b> compilation
 *                  flag was not used.
 * 
 * @note
 *       -# Though T is strictly TDB, it is usually more convenient to use
 *          TT, which makes no significant difference.  Julian centuries since
 *          J2000 is (JD - 2451545.0)/36525.
 *       -# The expression used is as adopted in IERS Conventions (2010) and
 *          is from Simon et al. (1994).  Arguments are in radians.
 *       -# L in this instance is the Mean Longitude of the Moon. OM is the 
 *          Mean longitude of the ascending node of the Moon.
 *       -# Status: Canonical model
 * 
 *  Test case:
 *     given input: T = 0.07995893223819302 Julian centuries since J2000
 *                  (MJD = 54465)
 *     expected output:  L = 2.291187512612069099 radians
 *                       LP = 6.212931111003726414 radians
 *                       F = 3.658025792050572989 radians
 *                       D = 4.554139562402433228 radians
 *                       OM = -0.5167379217231804489 radians
 * 
 * @version 2010 February 25
 * 
 * @cite iers2010, @cite simon94
 * 
 */
int iers2010::fundarg (const double& t,double& l,double& lp,double& f,double& d,double& om)
{

  // Set constants
  #ifdef USE_EXTERNAL_CONSTS
    /*constexpr double DAS2R   (DAS2R);     // Arcseconds to radians*/
    /*constexpr double TURNAS  (TURNAS);    // Arcseconds in a full circle*/
    constexpr double TWOPI   (D2PI);
  #else
    constexpr double DAS2R   ( 4.848136811095359935899141e-6 ); // Arcseconds to radians
    constexpr double TURNAS  ( 1296000e0 );                     // Arcseconds in a full circle
    constexpr double TWOPI   ( 6.283185307179586476925287e0 );  // 2*pi
  #endif
    
  #ifdef QUICK_EXIT
    static double t_previous  = .0e0;
    static double l_previous  = -999.0;
    static double lp_previous = -999.0;
    static double f_previous  = -999.0;
    static double d_previous  = -999.0;
    static double om_previous = -999.0;
    if ( fabs( t_previous-t ) < DATE_MAX_DIFF ) {
      l  = l_previous;
      lp = lp_previous;
      f  = f_previous;
      d  = d_previous;
      om = om_previous;
      return 1;
    }
  #endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  //  Compute the fundamental argument L.
  l = fmod (       485868.249036e0 +
           t*( 1717915923.2178e0 +
           t*(         31.8792e0 +
           t*(          0.051635e0 +
           t*(        - 0.00024470e0 )))), TURNAS ) * DAS2R;

  // Compute the fundamental argument LP.
  lp = fmod (       1287104.79305e0 +
              t*( 129596581.0481e0 +
              t*(       - 0.5532e0 +
              t*(         0.000136e0 +
              t*(       - 0.00001149e0 )))), TURNAS ) * DAS2R;

  // Compute the fundamental argument F.
  f  = fmod (       335779.526232e0 +
            t*( 1739527262.8478e0 +
            t*(       - 12.7512e0 +
            t*(       -  0.001037e0 +
            t*(          0.00000417e0 )))), TURNAS ) * DAS2R;

  // Compute the fundamental argument D.
  d = fmod (        1072260.70369e0 +
             t*( 1602961601.2090e0 +
             t*(        - 6.3706e0 +
             t*(          0.006593e0 +
             t*(        - 0.00003169e0 )))), TURNAS ) * DAS2R;

  // Compute the fundamental argument OM.
  om = fmod (       450160.398036e0 +
             t*( - 6962890.5431e0 +
             t*(         7.4722e0 +
             t*(         0.007702e0 +
             t*(       - 0.00005939e0 )))), TURNAS ) * DAS2R;

  #ifdef QUICK_EXIT
    t_previous  = t;
    l_previous  = l;
    lp_previous = lp;
    f_previous  = f;
    d_previous  = d;
    om_previous = om;
  #endif

  // Finished.
  return 0;
}