/**
 * @details This function computes the station tidal displacement
 *          caused by lunar and solar gravitational attraction (see References). 
 *          The computations are calculated by the following steps:<br>
 *          <b>Step 1):</b> General degree 2 and degree 3 corrections + CALL ST1IDIU 
 *          + CALL ST1ISEM + CALL ST1L1.<br>  
 *          <b>Step 2):</b> CALL STEP2DIU + CALL STEP2LON<br>
 *          It has been decided that the <b>Step 3</b> non-correction for permanent tide
 *          would not be applied in order to avoid a jump in the reference frame.
 *          This Step 3 must be added in order to get the non-tidal station position
 *          and to conform with the IAG Resolution.
 *          This function is a translation/wrapper for the fortran DEHANTTIDEINEL
 *          subroutine, found here : http://maia.usno.navy.mil/conv2010/software.html
 * 
*  Given:
*     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
*     XSUN          d(3)   Geocentric position of the Sun (Note 2)
*     XMON          d(3)   Geocentric position of the Moon (Note 2)
*     YR            i      Year (Note 3)
*     MONTH         i      Month (Note 3)
*     DAY           i      Day of Month (Note 3)
*     FHR           d      Hour in the day (Notes 3 and 4)
*
*  Returned:
*     DXTIDE        d(3)   Displacement vector (Note 5)
 * @param[in]  xsta    Geocentric position of the IGS station (Note 1)
 * @param[in]  fhr     Fractional hours in the day (Note 2)
 * @param[in]  t       Centuries since J2000
 * @param[out] xcorsta In phase and out of phase station corrections
 *                     for diurnal band (Note 4) 
 * 
 * @note
 *       -# The IGS station is in ITRF co-rotating frame. All coordinates are
 *          expressed in meters, as arrays, i.e. [x,y,z].
 *       -# The fractional hours in the day is computed as the hour + minutes/60.0
 *          + sec/3600.0.  The unit is expressed in Universal Time (UT).
 *       -# ----
 *       -# All coordinates are expressed in meters.
 *       -# Status: Class 1
 *       -# This fucnction is part of the package dehanttideinel, see
 *          ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/ 
 * 
 *  Test case:
 *     given input: XSTA(1) = 4075578.385D0 meters
 *                  XSTA(2) =  931852.890D0 meters
 *                  XSTA(3) = 4801570.154D0 meters 
 *                  FHR     = 0.00D0 hours
 *                  T       = 0.1059411362080767D0 Julian centuries
 *                  
 *     expected output:  XCORSTA(1) = 0.4193085327321284701D-02 meters
 *                       XCORSTA(2) = 0.1456681241014607395D-02 meters
 *                       XCORSTA(3) = 0.5123366597450316508D-02 meters
 * 
 * @version 2010 October  20
 * 
 * @cite iers2010,
 *       Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, "Tidal station
 *       displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477,
 * 
 */

#include <numeric>
#include <cmath>
#ifdef USE_EXTERNAL_CONSTS
  #include "gencon.hpp"
#endif

/**
 * @brief function to compute the scalar product of two vectors and their norms.
 * @param[in] x   vector of dimension 3
 * @param[in] y   vector of dimension 3
 * @param[out] r1 (euclidean) norm of vector x
 * @param[out] r2 (euclidean) norm of vector y
 * @return        scalar product of vectors x and y
 */
inline double _sprod_ (const double* x,const double* y,double& r1,double& r2) {
  r1 = ::sqrt ( std::inner_product (x,x+3,x,.0e0) );
  r2 = ::sqrt ( std::inner_product (y,y+3,y,.0e0) );
  return std::inner_product (x,x+3,y,.0e0);
}

int iers2010::dehanttideinel ()
{

  xcorsta[0] = xcorsta[1] = xcorsta[2] = .0e0;

  // Set constants
  #ifdef USE_EXTERNAL_CONSTS
    constexpr double PI   (DPI);                              // pi
  #else
    constexpr double PI   ( 3.1415926535897932384626433e0 );  // pi
  #endif

  /*----------------------------------------------------------------------  
   * NOMINAL SECOND DEGREE AND THIRD DEGREE LOVE NUMBERS AND SHIDA NUMBERS  
   *----------------------------------------------------------------------*/
  constexpr double h20 (0.6078e0), l20 (0.0847e0), h3 (0.292e0), l3 (0.015e0);

  /*----------------------------------------------------------------------  
   * SCALAR PRODUCT OF STATION VECTOR WITH SUN/MOON VECTOR  
   *----------------------------------------------------------------------*/
  double rsta,rsun,rmon;
  double scs = _sprod_ (xsta,xsun,rsta,rsun);
  double scm = _sprod_ (xsta,xmon,rsta,rmon); 
  double scsun = scs / rsta / rsun; 
  double scmon = scm / rsta / rmon;

  /*----------------------------------------------------------------------   
   * COMPUTATION OF NEW H2 AND L2  
   *----------------------------------------------------------------------*/
  double cosphi = ::sqrt ( xsta[0]*xsta[0]+xsta[1]*xsta[1] ) / rsta;
  double h2     = h20 - 0.0006d0*(1e0-3e0/2e0*cosphi*cosphi);
  double l2     = l20 + 0.0002e0*(1e0-3e0/2e0*cosphi*cosphi);

  // P2 term  
  double p2sun = 3e0*(h2/2e0-l2)*scsun*scsun-h2/2e0;
  double p2mon = 3e0*(h2/2e0-l2)*scmon*scmon-h2/2e0;

  // P3 term  
  double p3sun = 5e0/2e0*(h3-3e0*l3)*pow(scsun,3)+3e0/2e0*(l3-h3)*scsun;
  double p3mon = 5e0/2e0*(h3-3e0*l3)*pow(scmon,3)+3e0/2e0*(l3-h3)*scmon;

  /*----------------------------------------------------------------------  
   * TERM IN DIRECTION OF SUN/MOON VECTOR  
   *----------------------------------------------------------------------*/
  double x2sun = 3e0*l2*scsun;
  double x2mon = 3e0*l2*scmon;
  double x3sun = 3e0*l3/2e0*(5e0*scsun*scsun-1e0);
  double x3mon = 3e0*l3/2e0*(5e0*scmon*scmon-1e0);

  /*----------------------------------------------------------------------  
   * FACTORS FOR SUN/MOON USING IAU CURRENT BEST ESTIMATES (SEE REFERENCES) 
   *----------------------------------------------------------------------*/
  const double mass_ratio_sun  = 332946.0482e0;
  const double mass_ratio_moon = 0.0123000371e0;
  const double re              = 6378136.6e0;
  const double fac2sun         = mass_ratio_sun*re*pow(re/rsun,3);
  const double fac2mon         = mass_ratio_moon*re*pow(re/rmon,3);
  const double fac3sun         = fac2sun*(re/rsun);
  const double fac3mon         = fac2mon*(re/rmon);
  
  // TOTAL DISPLACEMENT
  double dxtide[] = {.0e0,.0e0,.0e0};
  for (int i=0;i<3;i++) 
    dxtide[i] = fac2sun*( x2sun*xsun[i]/rsun + p2sun*xsta[i]/rsta ) +
                fac2mon*( x2mon*xmon[i]/rmon + p2mon*xsta[i]/rsta ) +  
                fac3sun*( x3sun*xsun[i]/rsun + p3sun*xsta[i]/rsta ) +   
                fac3mon*( x3mon*xmon[i]/rmon + p3mon*xsta[i]/rsta );

  /*+---------------------------------------------------------------------  
   * CORRECTIONS FOR THE OUT-OF-PHASE PART OF LOVE NUMBERS (PART H_2^(0)I  
   * AND L_2^(0)I )  
   *----------------------------------------------------------------------*/

  // FIRST, FOR THE DIURNAL BAND
  st1idiu ( xsta,xsun,xmon,fac2sun,fac2mon,xcorsta );
  for (int i=0;i<3;i++) dxtide[i] += xcorsta[i];
  
  // SECOND, FOR THE SEMI-DIURNAL BAND       
  st1isem ( xsta,xsun,xmon,fac2sun,fac2mon,xcorsta );
  for (int i=0;i<3;i++) dxtide[i] += xcorsta[i];  


  /*+---------------------------------------------------------------------
   * CORRECTIONS FOR THE LATITUDE DEPENDENCE OF LOVE NUMBERS (PART L^(1) )  
   *----------------------------------------------------------------------*/
  st1l1 ( xsta,xsun,xmon,fac2sun,fac2mon,xcorsta );
  for (int i=0;i<3;i++) dxtide[i] += xcorsta[i];    
  
  // CONSIDER CORRECTIONS FOR STEP 2

  /*+---------------------------------------------------------------------  
   * CORRECTIONS FOR THE DIURNAL BAND:  
   * 
   *  FIRST, WE NEED TO KNOW THE DATE CONVERTED IN JULIAN CENTURIES 
   *        
   *   1) CALL THE SUBROUTINE COMPUTING THE JULIAN DATE 
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 
   ///      CALL CAL2JD ( YR, MONTH, DAY, JJM0, JJM1, STATUT )
   ///      FHRD = FHR/24.D0
  /*    17 May 2013 Corrected bug as noted in header                     */
   ///      T=((JJM0-2451545.0D0)+JJM1+FHRD)/36525D0
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   *   2) CALL THE SUBROUTINE COMPUTING THE CORRECTION OF UTC TIME  
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
   ///     CALL DAT ( YR, MONTH, DAY, FHRD, DTT, STATUT )
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
   ///      DTT = DTT + 32.184D0
  /*    CONVERSION OF T IN TT TIME                                        */
   ///     T=T+DTT/(3600.0D0*24.0D0*36525D0)

  //  SECOND, WE CAN CALL THE SUBROUTINE STEP2DIU, FOR THE DIURNAL BAND
  //+ CORRECTIONS, (in-phase and out-of-phase frequency dependence):
  step2diu ( xsta,fhr,t,xcorsta );
  for (int i=0;i<3;i++) dxtide[i] += xcorsta[i]; 
  
  //  CORRECTIONS FOR THE LONG-PERIOD BAND,
  //+ (in-phase and out-of-phase frequency dependence):  
  step2lon ( xsta,t,xcorsta );
  for (int i=0;i<3;i++) dxtide[i] += xcorsta[i];  
    
  // CONSIDER CORRECTIONS FOR STEP 3

  /*----------------------------------------------------------------------
   * UNCORRECT FOR THE PERMANENT TIDE  
   *-----------------------------------------------------------------------*/
/*
  double sinphi = xsta[2]/rsta;
  double cosphi = sqrt( xsta[0]*xsta[0]+xsta[1]*xsta[1] ) / rsta;
  double cosla  = xsta[0]/cosphi/rsta;
  double sinla  = xsta[1]/cosphi/rsta;
  double dr     = -sqrt(5e0/4e0/PI)*h2*0.31460e0*(3e0/2e0*sinphi*sinphi-0.5e0);
  double dn     = -sqrt(5e0/4e0/PI)*l2*0.31460e0*3e0*cosphi*sinphi;
  dxtide[0]    += -dr*cosla*cosphi+dn*cosla*sinphi;
  dxtide[1]    += -dr*sinla*cosphi+dn*sinla*sinphi;
  dxtide[2]    += -dr*sinphi      -dn*cosphi;
*/

  // Finished.
  return 0;
}
