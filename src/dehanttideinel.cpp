#include "iers2010.hpp"
#include "dehanttideinel.hpp"

/**
 * @details This function computes the station tidal displacement
 *          caused by lunar and solar gravitational attraction (see References).
 *          The computations are calculated by the following steps:<br>
 *          <b>Step 1):</b> General degree 2 and degree 3 corrections 
 *          + CALL ST1IDIU + CALL ST1ISEM + CALL ST1L1.<br>  
 *          <b>Step 2):</b> CALL STEP2DIU + CALL STEP2LON<br>
 *          It has been decided that the <b>Step 3</b> non-correction for 
 *          permanent tide would not be applied in order to avoid a jump in the
 *          reference frame. This Step 3 must be added in order to get the 
 *          non-tidal station position and to conform with the IAG Resolution.
 *          This function is a translation/wrapper for the fortran 
 *          DEHANTTIDEINEL subroutine, found here : 
 *          http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  xsta   Geocentric position of the IGS station (Note 1)
 * @param[in]  xsun   Geocentric position of the Sun (Note 2)
 * @param[in]  xmon   Geocentric position of the Moon (Note 2)
 * @param[in]  yr     Year (Note 3)
 * @param[in]  month  Month (Note 3)
 * @param[in]  day    Day of Month (Note 3)
 * @param[in]  fhr    Hour in the day (Notes 3 and 4)
 * @param[out] dxtide Displacement vector (Note 5) 
 * 
 * @note
 *     -# The IGS station is in ITRF co-rotating frame.  All coordinates,
 *        X, Y, and Z, are expressed in meters. 
 *     -# The position is in Earth Centered Earth Fixed (ECEF) frame.  All
 *        coordinates are expressed in meters.
 *     -# The values are expressed in Coordinated Universal Time (UTC).
 *     -# The fractional hours in the day is computed as the hour + minutes/60.0
 *        + sec/3600.0.
 *     -# The displacement vector is in the geocentric ITRF.  All components are
 *        expressed in meters.
 *     -# Parameters jc1 and jc2 constitute the date as Julian Centuries in TT 
 *        time scale. The actual date is given by the addition jc1+jc2. 
 *        Either jc1 or jc2 can be set to zero.
 *     -# Status: Class 1
 *     -# This fucnction is part of the package dehanttideinel, see
 *        ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/ 
 * 
 * @verbatim
 *  Test case:
 *     given input: XSTA(1) = 4075578.385D0 meters
 *                  XSTA(2) =  931852.890D0 meters
 *                  XSTA(3) = 4801570.154D0 meters   
 *                  XSUN(1) = 137859926952.015D0 meters
 *                  XSUN(2) = 54228127881.4350D0 meters
 *                  XSUN(3) = 23509422341.6960D0 meters
 *                  XMON(1) = -179996231.920342D0 meters
 *                  XMON(2) = -312468450.131567D0 meters
 *                  XMON(3) = -169288918.592160D0 meters
 *                  YR      = 2009
 *                  MONTH   = 4
 *                  DAY     = 13
 *                  FHR     = 0.00D0 hour
 *                  
 *     expected output:  DXTIDE(1) = 0.7700420357108125891D-01 meters
 *                       DXTIDE(2) = 0.6304056321824967613D-01 meters
 *                       DXTIDE(3) = 0.5516568152597246810D-01 meters
 *
 *  Test case:
 *     given input: XSTA(1) =  1112189.660D0 meters
 *                  XSTA(2) = -4842955.026D0 meters
 *                  XSTA(3) =  3985352.284D0 meters   
 *                  XSUN(1) = -54537460436.2357D0 meters
 *                  XSUN(2) =  130244288385.279D0 meters
 *                  XSUN(3) =  56463429031.5996D0 meters
 *                  XMON(1) =  300396716.912D0 meters
 *                  XMON(2) =  243238281.451D0 meters
 *                  XMON(3) =  120548075.939D0 meters
 *                  YR      = 2012
 *                  MONTH   = 7
 *                  DAY     = 13
 *                  FHR     = 0.00D0 hour
 *                  
 *     expected output:  DXTIDE(1) = -0.2036831479592075833D-01 meters
 *                       DXTIDE(2) =  0.5658254776225972449D-01 meters
 *                       DXTIDE(3) = -0.7597679676871742227D-01 meters
 * @endverbatim
 * 
 * @version 2013 May      17
 *
 * @cite iers2010,
 *
 *     Groten, E., 2000, Geodesists Handbook 2000, Part 4,
 *     http://www.gfy.ku.dk/~iag/HB2000/part4/groten.htm. See also
 *     ''Parameters of Common Relevance of Astronomy, Geodesy, and
 *     Geodynamics," J. Geod., 74, pp. 134-140
 *
 *     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
 *     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
 *
 *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
 *     IERS Technical Note No. 36, BKG (2010)
 *     
 *     Pitjeva, E. and Standish, E. M., 2009, ''Proposals for the masses
 *     of the three largest asteroids, the Moon-Earth mass ratio and the
 *     Astronomical Unit," Celest. Mech. Dyn. Astr., 103, pp. 365-372
 *
 *     Ries, J. C., Eanes, R. J., Shum, C. K. and Watkins, M. M., 1992,
 *     ''Progress in the Determination of the Gravitational Coefficient
 *     of the Earth," Geophys. Res. Lett., 19(6), pp. 529-531
 * 
 */
int iers2010::dehanttideinel (const double* xsta,const double* xsun, 
    const double* xmon,const int& yr,const int& month, 
    const int& day,const double& fhr,double* dxtide)
{
    
    double xcorsta[] = {.0e0,.0e0,.0e0};
    
    /*----------------------------------------------------------------------  
    * NOMINAL SECOND DEGREE AND THIRD DEGREE LOVE NUMBERS AND SHIDA NUMBERS  
    *----------------------------------------------------------------------*/
    constexpr double h20 (0.6078e0), l20 (0.0847e0), 
              h3 (0.292e0), l3 (0.015e0);
    
    /*----------------------------------------------------------------------  
    * SCALAR PRODUCT OF STATION VECTOR WITH SUN/MOON VECTOR  
    *----------------------------------------------------------------------*/
    double rsta,rsun,rmon;
    double scs   ( iers2010::dtel::_sprod_ (xsta,xsun,rsta,rsun) );
    double scm   ( iers2010::dtel::_sprod_ (xsta,xmon,rsta,rmon) );
    double scsun ( scs / rsta / rsun );
    double scmon ( scm / rsta / rmon );
    
    /*----------------------------------------------------------------------   
    * COMPUTATION OF NEW H2 AND L2  
    *----------------------------------------------------------------------*/
    double cosphi ( ::sqrt ( xsta[0]*xsta[0]+xsta[1]*xsta[1] ) / rsta );
    double h2     ( h20 - 0.0006e0*(1e0-3e0/2e0*cosphi*cosphi) );
    double l2     ( l20 + 0.0002e0*(1e0-3e0/2e0*cosphi*cosphi) );
    
    // P2 term  
    double p2sun ( 3e0*(h2/2e0-l2)*scsun*scsun-h2/2e0 );
    double p2mon ( 3e0*(h2/2e0-l2)*scmon*scmon-h2/2e0 );
    
    // P3 term  
    double p3sun ( 5e0/2e0*(h3-3e0*l3)*pow(scsun,3)+3e0/2e0*(l3-h3)*scsun );
    double p3mon ( 5e0/2e0*(h3-3e0*l3)*pow(scmon,3)+3e0/2e0*(l3-h3)*scmon );
    
    /*----------------------------------------------------------------------  
    * TERM IN DIRECTION OF SUN/MOON VECTOR  
    *----------------------------------------------------------------------*/
    double x2sun ( 3e0*l2*scsun );
    double x2mon ( 3e0*l2*scmon );
    double x3sun ( 3e0*l3/2e0*(5e0*scsun*scsun-1e0) );
    double x3mon ( 3e0*l3/2e0*(5e0*scmon*scmon-1e0) );
    
    /*----------------------------------------------------------------------  
    * FACTORS FOR SUN/MOON USING IAU CURRENT BEST ESTIMATES (SEE REFERENCES) 
    *----------------------------------------------------------------------*/
    const double mass_ratio_sun  ( 332946.0482e0 );
    const double mass_ratio_moon ( 0.0123000371e0 );
    const double re              ( 6378136.6e0 );
    const double fac2sun         ( mass_ratio_sun*re*pow(re/rsun,3) );
    const double fac2mon         ( mass_ratio_moon*re*pow(re/rmon,3) );
    const double fac3sun         ( fac2sun*(re/rsun) );
    const double fac3mon         ( fac2mon*(re/rmon) );
    
    // TOTAL DISPLACEMENT
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
    iers2010::dtel::st1idiu (xsta,xsun,xmon,fac2sun,fac2mon,xcorsta);
    for (int i=0;i<3;i++)
        dxtide[i] += xcorsta[i];
    
    // SECOND, FOR THE SEMI-DIURNAL BAND       
    iers2010::dtel::st1isem (xsta,xsun,xmon,fac2sun,fac2mon,xcorsta);
    for (int i=0;i<3;i++)
        dxtide[i] += xcorsta[i];  
    
    
    /*+---------------------------------------------------------------------
    * CORRECTIONS FOR THE LATITUDE DEPENDENCE OF LOVE NUMBERS (PART L^(1) )  
    *----------------------------------------------------------------------*/
    iers2010::dtel::st1l1 (xsta,xsun,xmon,fac2sun,fac2mon,xcorsta);
    for (int i=0;i<3;i++)
        dxtide[i] += xcorsta[i];
    
    // CONSIDER CORRECTIONS FOR STEP 2
    
    /*+---------------------------------------------------------------------  
    * CORRECTIONS FOR THE DIURNAL BAND:  
    * 
    *  FIRST, WE NEED TO KNOW THE DATE CONVERTED IN JULIAN CENTURIES 
    *        
    *   1) CALL THE SUBROUTINE COMPUTING THE JULIAN DATE 
    *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    double jjm0,jjm1;
    iers2010::dtel::cal2jd (yr,month,day,jjm0,jjm1);
    double fhrd ( fhr / 24.0e0 );
    /*    17 May 2013 Corrected bug as noted in header                     */
    double t ( ((jjm0-2451545.0e0)+jjm1+fhrd)/36525e0 );
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    *   2) CALL THE SUBROUTINE COMPUTING THE CORRECTION OF UTC TIME  
    *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    double dtt;
    iers2010::dtel::dat (yr,month,day,fhrd,dtt);
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    dtt += 32.184e0;
    /*    CONVERSION OF T IN TT TIME                                        */
    t += dtt/(3600.0e0*24.0e0*36525e0);
    
    //  SECOND, WE CAN CALL THE SUBROUTINE STEP2DIU, FOR THE DIURNAL BAND
    //+ CORRECTIONS, (in-phase and out-of-phase frequency dependence):
    iers2010::dtel::step2diu (xsta,fhr,t,xcorsta);
    for (int i=0;i<3;i++)
        dxtide[i] += xcorsta[i];
    
    //  CORRECTIONS FOR THE LONG-PERIOD BAND,
    //+ (in-phase and out-of-phase frequency dependence):  
    iers2010::dtel::step2lon (xsta,t,xcorsta);
    for (int i=0;i<3;i++)
        dxtide[i] += xcorsta[i];
    
    // CONSIDER CORRECTIONS FOR STEP 3
    
    /*----------------------------------------------------------------------
    * UNCORRECT FOR THE PERMANENT TIDE  
    *-----------------------------------------------------------------------*/
    /*
    * double sinphi = xsta[2]/rsta;
    * double cosphi = sqrt( xsta[0]*xsta[0]+xsta[1]*xsta[1] ) / rsta;
    * double cosla  = xsta[0]/cosphi/rsta;
    * double sinla  = xsta[1]/cosphi/rsta;
    * double dr     = -sqrt(5e0/4e0/PI)*h2*0.31460e0*(3e0/2e0*sinphi*sinphi-0.5e0);
    * double dn     = -sqrt(5e0/4e0/PI)*l2*0.31460e0*3e0*cosphi*sinphi;
    * dxtide[0]    += -dr*cosla*cosphi+dn*cosla*sinphi;
    * dxtide[1]    += -dr*sinla*cosphi+dn*sinla*sinphi;
    * dxtide[2]    += -dr*sinphi      -dn*cosphi;
    */
    
    // Finished.
    return 0;
}