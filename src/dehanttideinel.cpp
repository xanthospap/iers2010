#include "iers2010.hpp"

inline double
sprod(const double* x1, const double* x2, double& r1, double& r2)
{
    r1 = std::sqrt(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2]);
    r2 = std::sqrt(x2[0]*x2[0]+x2[1]*x2[1]+x2[2]*x2[2]);
    return x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
}

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
 * @param[in]  xsta   Geocentric position of the station (Note 1)
 * @param[in]  xsun   Geocentric position of the Sun (Note 2)
 * @param[in]  xmon   Geocentric position of the Moon (Note 2)
 * @param[in]  yr     Year (Note 3)
 * @param[in]  month  Month (Note 3)
 * @param[in]  day    Day of Month (Note 3)
 * @param[in]  fhr    Hour in the day (Notes 3 and 4)
 * @param[out] dxtide Displacement vector (Note 5)
 * @return            Always 0.
 * 
 * @note
 *     -# The station is in ITRF co-rotating frame.  All coordinates,
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
 * @version 19.12.2016
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

int
iers2010::dehanttideinel(const double* xsta,const double* xsun, 
    const double* xmon, int yr, int month, int day, double fhr,
    double* dxtide)
{
    
    double xcorsta[] = {.0e0,.0e0,.0e0};
    
    /*----------------------------------------------------------------------  
     * nominal second degree and third degree love numbers and shida numbers  
     *----------------------------------------------------------------------*/
    constexpr double h20 { 0.6078e0 },
                     l20 { 0.0847e0 }, 
                     h3  { 0.292e0  },
                     l3  { 0.015e0  };
    
    /*----------------------------------------------------------------------  
     * scalar product of station vector with sun/moon vector  
     *----------------------------------------------------------------------*/
    double rsta, rsun, rmon;
    double scs   { sprod(xsta, xsun, rsta, rsun) };
    double scm   { sprod(xsta, xmon, rsta, rmon) };
    double scsun { scs / rsta / rsun };
    double scmon { scm / rsta / rmon };
    
    /*----------------------------------------------------------------------   
     * computation of new h2 and l2  
     *----------------------------------------------------------------------*/
    double cosphi { std::sqrt(xsta[0]*xsta[0]+xsta[1]*xsta[1] ) / rsta };
    double h2     { h20 - 0.0006e0*(1e0-3e0/2e0*cosphi*cosphi) };
    double l2     { l20 + 0.0002e0*(1e0-3e0/2e0*cosphi*cosphi) };
    
    // P2 term  
    double p2sun { 3e0*(h2/2e0-l2)*scsun*scsun-h2/2e0 };
    double p2mon { 3e0*(h2/2e0-l2)*scmon*scmon-h2/2e0 };
    
    // P3 term  
    double p3sun { 5e0/2e0*(h3-3e0*l3)*pow(scsun,3)+3e0/2e0*(l3-h3)*scsun };
    double p3mon { 5e0/2e0*(h3-3e0*l3)*pow(scmon,3)+3e0/2e0*(l3-h3)*scmon };
    
    /*----------------------------------------------------------------------  
     * term in direction of sun/moon vector  
     *----------------------------------------------------------------------*/
    double x2sun { 3e0*l2*scsun };
    double x2mon { 3e0*l2*scmon };
    double x3sun { 3e0*l3/2e0*(5e0*scsun*scsun-1e0) };
    double x3mon { 3e0*l3/2e0*(5e0*scmon*scmon-1e0) };
    
    /*----------------------------------------------------------------------  
     * factors for sun/moon using iau current best estimates (see references) 
     *----------------------------------------------------------------------*/
    const double mass_ratio_sun  { 332946.0482e0  };
    const double mass_ratio_moon { 0.0123000371e0 };
    const double re              { 6378136.6e0    };
    const double fac2sun         { mass_ratio_sun*re*pow(re/rsun,3) };
    const double fac2mon         { mass_ratio_moon*re*pow(re/rmon,3) };
    const double fac3sun         { fac2sun*(re/rsun) };
    const double fac3mon         { fac2mon*(re/rmon) };
    
    // total displacement
    for (int i=0; i<3; i++) {
        dxtide[i] = fac2sun*( x2sun*xsun[i]/rsun + p2sun*xsta[i]/rsta ) +
                    fac2mon*( x2mon*xmon[i]/rmon + p2mon*xsta[i]/rsta ) +
                    fac3sun*( x3sun*xsun[i]/rsun + p3sun*xsta[i]/rsta ) +
                    fac3mon*( x3mon*xmon[i]/rmon + p3mon*xsta[i]/rsta );
    }
    
    /*+---------------------------------------------------------------------  
    * corrections for the out-of-phase part of love numbers (part h_2^(0)i  
    * and l_2^(0)i )  
    *----------------------------------------------------------------------*/
    
    // first, for the diurnal band
    iers2010::dhtide::st1idiu(xsta, xsun, xmon, fac2sun, fac2mon, xcorsta);
    dxtide[0] += xcorsta[0];
    dxtide[1] += xcorsta[1];
    dxtide[2] += xcorsta[2];

    // second, for the semi-diurnal band       
    iers2010::dhtide::st1isem(xsta, xsun, xmon, fac2sun, fac2mon, xcorsta);
    dxtide[0] += xcorsta[0];
    dxtide[1] += xcorsta[1];
    dxtide[2] += xcorsta[2];
    
    /*---------------------------------------------------------------------
     * corrections for the latitude dependence of love numbers (part l^(1) )  
     *----------------------------------------------------------------------*/
    iers2010::dhtide::st1l1(xsta, xsun, xmon, fac2sun, fac2mon, xcorsta);
    dxtide[0] += xcorsta[0];
    dxtide[1] += xcorsta[1];
    dxtide[2] += xcorsta[2];
    
    // CONSIDER CORRECTIONS FOR STEP 2
    
    /*---------------------------------------------------------------------  
     * corrections for the diurnal band:  
     * 
     *  first, we need to know the date converted in julian centuries 
     *        
     *   1) call the subroutine computing the julian date 
     *---------------------------------------------------------------------*/
    double jjm0, jjm1;
    iers2010::dhtide::cal2jd(yr, month, day, jjm0, jjm1);
    double fhrd { fhr / 24.0e0 };
    // 17 May 2013 Corrected bug as noted in header
    double t { ((jjm0-2451545.0e0)+jjm1+fhrd)/36525e0 };

    // 
    //  2) call the subroutine computing the correction of utc time  
    //
    double dtt;
    iers2010::dhtide::dat(yr, month, day, fhrd, dtt);
    dtt += 32.184e0;
    // conversion of t in tt time
    t += dtt/(3600.0e0*24.0e0*36525e0);
    
    //  second, we can call the subroutine step2diu, for the diurnal band
    //+ corrections, (in-phase and out-of-phase frequency dependence):
    iers2010::dhtide::step2diu(xsta, fhr, t, xcorsta);
    dxtide[0] += xcorsta[0];
    dxtide[1] += xcorsta[1];
    dxtide[2] += xcorsta[2];
    
    //  corrections for the long-period band,
    //+ (in-phase and out-of-phase frequency dependence):  
    iers2010::dhtide::step2lon(xsta, t, xcorsta);
    dxtide[0] += xcorsta[0];
    dxtide[1] += xcorsta[1];
    dxtide[2] += xcorsta[2];
    
    // CONSIDER CORRECTIONS FOR STEP 3
    
    /*----------------------------------------------------------------------
     * uncorrect for the permanent tide  
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
