/**
 * @details This subroutine computes the effects of the free core nutation.  
 *          Please note that the table is updated each year (see Note 4).
 *          The parameter N needs to be incremented for each additional
 *          year in the table.
 * 
 * @param[in] mjd   Modified Julian Date, TDB (Note 1)
 * @param[out]  x   CIP offset x component, in microas (Note 2)
 * @param[out]  y   CIP offset y component, in microas (Note 2)
 * @param[out] dx   Uncertainty of x component, in microas (Note 3)
 * @param[out] dy   Uncertainty of y component, in microas (Note 3) 
 * @return          An integer value which can be either 0 (to denote
 *                  that a new value for pm has been computed) or 1 
 *                  (to denote quick return). The returned integer has
 *                  nothing to do with the status of the function and
 *                  has no meaning if the <b>QUICK_EXIT</b> compilation
 *                  flag was not used.
 * 
 * @warning Contains data table to be updated each year, see http://syrte.obspm.fr/~lambert/fcn/
 * 
 * @note 
 *      -#  Though the Modified Julian Date (MJD) is strictly TDB, it is
 *          usually more convenient to use TT, which makes no significant
 *          difference.
 *      -#  CIP is the Celestial Intermediate Pole.  The expression
 *          used is given in microarcseconds.
 *      -#  The expression used is given in microarcseconds.
 *      -#  The updated table is maintained at the website
 *          http://syrte.obspm.fr/~lambert/fcn/.
 *      -#  Status: Class 3 model
 * 
 *  Test case: (NOT UPDATED FOR 2013 TABLE, TO BE DONE)
 *     given input: MJD = 54790D0   Modified Julian Date, TDB
 *                  
 *     expected output:  X = -176.8012290066270680D0 microarcseconds
 *                       Y = -93.51855308903756736D0 microarcseconds
 *                       dX = 3.745573770491803067D0 microarcseconds
 *                       dY = 3.745573770491803067D0 microarcseconds
 * 
 * @version 2013 December 19
 * 
 * @cite iers2010
 * 
 */

#include "iers2010.hpp"
#ifdef USE_EXTERNAL_CONSTS
  #include "gencon.hpp"
#endif

int iers2010::fcnnut (const double& mjd,double& x,double& y,double& dx,double& dy) 
{

  // Set constants
  static const int N = 30;

  #ifdef USE_EXTERNAL_CONSTS
    constexpr double PI   (DPI);
  #else
    constexpr double PI   ( 3.14159265358979323846e0 );     /// pi
  #endif
  // Mean prediction error
  constexpr double mpe    ( 0.1325e0 );                     /// microarcseconds per day
  // FCN parameters
  constexpr double per    ( -430.21e0 );                    /// period in days
  double phi    ( (2e0*PI/per)*(mjd-51544.5e0) ); /// phase in radians

  // Block data of amplitudes for X (microas)
  static const struct {
    double date,mxc,mxs,msx;
    inline double xc () const { return  mxc; }
    inline double xs () const { return  mxs; }
    inline double sx () const { return  msx; }    
    inline double yc () const { return  mxs; }
    inline double ys () const { return -mxc; }
    inline double sy () const { return  msx; }
  } table [] = {
    //DATA ((TABLE(I,J),I=1,4),J=1,N) /
    {45700.e0,     4.55e0,   -36.58e0,    19.72e0}, ///! 1984.0
    {46066.e0,  -141.82e0,  -105.35e0,    11.12e0}, ///! 1985.0
    {46431.e0,  -246.56e0,  -170.21e0,     9.47e0}, ///! 1986.0
    {46796.e0,  -281.89e0,  -159.24e0,     8.65e0}, ///! 1987.0
    {47161.e0,  -255.05e0,   -43.58e0,     8.11e0}, ///! 1988.0
    {47527.e0,  -210.46e0,   -88.56e0,     7.31e0}, ///! 1989.0
    {47892.e0,  -187.79e0,   -57.35e0,     6.41e0}, ///! 1990.0
    {48257.e0,  -163.01e0,    26.26e0,     5.52e0}, ///! 1991.0
    {48622.e0,  -145.53e0,    44.65e0,     4.80e0}, ///! 1992.0
    {48988.e0,  -145.12e0,    51.49e0,     5.95e0}, ///! 1993.0
    {49353.e0,  -109.93e0,    16.87e0,     9.45e0}, ///! 1994.0
    {49718.e0,   -87.30e0,     5.36e0,     8.25e0}, ///! 1995.0
    {50083.e0,   -90.61e0,     1.52e0,     7.67e0}, ///! 1996.0
    {50449.e0,   -94.73e0,    35.35e0,     4.40e0}, ///! 1997.0
    {50814.e0,   -67.52e0,    27.57e0,     3.40e0}, ///! 1998.0
    {51179.e0,   -44.11e0,   -14.31e0,     3.45e0}, ///! 1999.0
    {51544.e0,     5.21e0,   -74.87e0,     3.26e0}, ///! 2000.0
    {51910.e0,    70.37e0,  -129.66e0,     2.86e0}, ///! 2001.0
    {52275.e0,    86.47e0,  -127.84e0,     2.75e0}, ///! 2002.0
    {52640.e0,   110.44e0,   -42.73e0,     2.59e0}, ///! 2003.0
    {53005.e0,   114.78e0,    -0.13e0,     2.53e0}, ///! 2004.0
    {53371.e0,   132.96e0,    -4.78e0,     2.72e0}, ///! 2005.0
    {53736.e0,   157.36e0,    28.63e0,     2.19e0}, ///! 2006.0
    {54101.e0,   160.40e0,    58.87e0,     1.87e0}, ///! 2007.0
    {54466.e0,   156.76e0,   101.24e0,     1.74e0}, ///! 2008.0
    {54832.e0,   142.99e0,   143.01e0,     1.89e0}, ///! 2009.0
    {55197.e0,    33.70e0,   184.46e0,     1.95e0}, ///! 2010.0
    {55562.e0,     0.76e0,   253.70e0,     1.14e0}, ///! 2011.0
    {55927.e0,    25.47e0,   271.66e0,     1.07e0}, ///! 2012.0
    {56293.e0,   113.42e0,   256.50e0,     1.86e0}  ///! 2013.0
  };

  // Prediction of the amplitude at the input date
  double axc,axs,ayc,ays;
  if ( mjd<=table[0].date ) {
    axc = table[0].xc();
    axs = table[0].xs();
    ayc = table[0].yc();
    ays = table[0].ys();
  } else if ( mjd>=table[N-1].date ) {
    axc = table[N-1].xc();
    axs = table[N-1].xs();
    ayc = table[N-1].yc();
    ays = table[N-1].ys();
  } else {
    for ( int i=0;i<N-1;i++ ) {
      if ( mjd>=table[i].date ) {
        double t = mjd - table[i].date;
        double dt   = table[i+1].date - table[i].date;
        double daxc = table[i+1].xc() - table[i].xc();
        double daxs = table[i+1].xs() - table[i].xs();
        double dayc = table[i+1].yc() - table[i].yc();
        double days = table[i+1].ys() - table[i].ys();
        axc  = table[i].xc() + ( daxc / dt )*t;
        axs  = table[i].xs() + ( daxs / dt )*t;
        ayc  = table[i].yc() + ( dayc / dt )*t;
        ays  = table[i].ys() + ( days / dt )*t;
      }
    }
  }

  // Computation of X and Y
  double cosphi ( cos(phi) ),sinphi ( sin(phi) ); 
  x = axc * cosphi - axs * sinphi;  /// microas
  y = ayc * cosphi - ays * sinphi;  /// microas

  // Prediction of the uncertainty at the input date
  if ( mjd<=table[0].date ) {
    axc = table[0].sx() + mpe * ( table[0].date - mjd );
    axs = table[0].sx() + mpe * ( table[0].date - mjd );
    ayc = table[0].sy() + mpe * ( table[0].date - mjd );
    ays = table[0].sy() + mpe * ( table[0].date - mjd );
  } else if ( mjd>=table[N-1].date ) {
    axc = table[N-1].sx() + mpe * ( mjd - table[N-1].date );
    axs = table[N-1].sx() + mpe * ( mjd - table[N-1].date );
    ayc = table[N-1].sy() + mpe * ( mjd - table[N-1].date );
    ays = table[N-1].sy() + mpe * ( mjd - table[N-1].date );
  } else {
    for ( int i=0;i<N-1;i++ ) {
      if ( mjd>=table[i].date ) {
        double t = mjd - table[i].date;
        double dt   = table[i+1].date - table[i].date;
        double daxc = table[i+1].sx() - table[i].sx();
        double daxs = table[i+1].sx() - table[i].sx();
        double dayc = table[i+1].sy() - table[i].sy();
        double days = table[i+1].sy() - table[i].sy();
        axc = fabs( table[i].sx() + ( daxc /dt )*t );
        axs = fabs( table[i].sx() + ( daxs /dt )*t );
        ayc = fabs( table[i].sy() + ( dayc /dt )*t );
        ays = fabs( table[i].sy() + ( days /dt )*t );
      }
    }
  }

  // Computation of the uncertainties
  dx = axc + axs; /// microas
  dy = ayc + ays; /// microas

  // Finished
  return 0;
}