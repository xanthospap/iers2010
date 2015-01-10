#include "hardisp.hpp"

/**
 * @details This subroutine returns the frequency and phase of a tidal
 *          constituent when its Doodson number is given as input. 
 * 
 * @param[in]  idood Doodson number of a tidal constituent
 * @param[in]  itm   Date as integer array in UTC. The format should be:
 *                   [year,day_of_year,hours,minutes,seconds].
 * @param[out] freq  Frequency of a tidal constituent
 * @param[out] phase Phase of a tidal constituent (Note 1)
 * 
 * @note
 *     -# The phases must be decreased by 90 degrees if the sum of the order 
 *        and the species number is odd (as for the 2nd degree diurnals, and 
 *        3rd degree low frequency and semidiurnals).
 *        These phases may need further adjustment to allow for the spherical
 *        harmonic normalization used; e.g. for that used for the potential
 *        by Cartwright and Tayler, 180 degrees must be added for (species,
 *        order) = (1,2), (1,3), or (3,3).
 *     -# Status:  Class 1 model
 * 
 * @todo Generate a version of this function to use Julian Dates. Date manipulation/
 *       transformation should not be used inside a sub-function.
 * 
 * @verbatim
 *   Test case:
 *     given input: For June 25, 2009 0 Hr 0 Min, M2 tide
 *                  DATA IDOOD = / 2, 0, 0, 0, 0, 0 /  
 *
 *     expected output: FREQ = 1.93227361605688D0
 *                      PHASE = 303.362015399442D0
 * @endverbatim
 * 
 * @cite iers2010
 * 
 * @version 2013 September 11
 * 
 */
int iers2010::hisp::dtfrph (const int* idood,const int* itm,double& freq,double& phase)
{
  static int last_date[] = {0,0,0,0,0};
  static double  d[6];
  static double dd[6];
  int itm2[6];
  
  /*------------------------------------------------------------------------
   *  Test to see if time has changed; if so, set the phases and frequencies
   *  for each of the Doodson arguments
   *------------------------------------------------------------------------*/
  bool initial = false;
  for (int i=0;i<5;i++)
    if (last_date[i] != itm[i])
      initial = true;
    
  if (initial) { /* Need to re-copute vectors d and dd based on date */
    
    // Convert times to Julian days (UT) then to Julian centuries from J2000.0 (ET)
    
    // This will only use itm[0:2] and itm2[0:3], ie. itm1[0]=year and itm[1]=day of year
    iers2010::hisp::toymd (itm,itm2);
    int jd = iers2010::hisp::juldat(itm2);
    double dayfr = itm[2]/24.0e0 + itm[3]/1440.0e0 + itm[4]/84600.0e0;
    double year  = itm[0]+(itm[1]+dayfr) / 
                  (365.0e0+(double)iers2010::hisp::leap(itm[0]));
    double delta = iers2010::hisp::etutc(year,delta);
    double djd   = (double)jd - 0.5e0 + dayfr;
    double t     = (djd - 2451545.0e0 + delta/86400.0e0)/36525.0e0;
    
    // IERS expressions for the Delaunay arguments, in degrees
    
    double f1 = 134.9634025100e0 +
         t*( 477198.8675605000e0 +
         t*(      0.0088553333e0 +
         t*(      0.0000143431e0 +
         t*(     -0.0000000680e0 ))));
    double f2 = 357.5291091806e0 +
         t*(  35999.0502911389e0 +
         t*(     -0.0001536667e0 +
         t*(      0.0000000378e0 +
         t*(     -0.0000000032e0 ))));
    double f3 = 93.2720906200e0 +
        t*( 483202.0174577222e0 +
        t*(     -0.0035420000e0 +
        t*(     -0.0000002881e0 +
        t*(      0.0000000012e0 ))));
    double f4 = 297.8501954694e0 +
         t*( 445267.1114469445e0 +
         t*(     -0.0017696111e0 +
         t*(      0.0000018314e0 +
         t*(     -0.0000000088e0 ))));
    double f5 = 125.0445550100e0 +
         t*(  -1934.1362619722e0 +
         t*(      0.0020756111e0 +
         t*(      0.0000021394e0 +
         t*(     -0.0000000165e0 ))));
         
    // Convert to Doodson (Darwin) variables
         
    d[0] = 360.0e0*dayfr - f4;
    d[1] = f3 + f5;
    d[2] = d[1] - f4;
    d[3] = d[1] - f1;
    d[4] = -f5;
    d[5] = d[2] - f2;
    
    //  Find frequencies of Delauney variables (in cycles/day), and from these
    //+ the same for the Doodson arguments
    
    double fd1 =  0.0362916471e0 + 0.0000000013e0*t;
    double fd2 =  0.0027377786e0;
    double fd3 =  0.0367481951e0 - 0.0000000005e0*t;
    double fd4 =  0.0338631920e0 - 0.0000000003e0*t;
    double fd5 = -0.0001470938e0 + 0.0000000003e0*t;
    dd[0] = 1.0e0 - fd4;
    dd[1] = fd3 + fd5;
    dd[2] = dd[1] - fd4;
    dd[3] = dd[1] - fd1;
    dd[4] = -fd5;
    dd[5] = dd[2] - fd2;
  }
  
  
  //  Compute phase and frequency of the given tidal constituent
  
  freq  = 0.0e0;
  phase = 0.0e0;
  for (int i=0;i<6;i++) {
    freq  += idood[i]*dd[i];
    phase += idood[i]*d[i];
  }
  
  // Adjust phases so that they fall in the positive range 0 to 360
  phase = fmod (phase,360.0e0);
  if (phase < 0.0e0) phase += 360.0e0;
  
  // copy the just used date to the static one for next use of function
  if (initial) std::copy (itm,itm+5,last_date);
  
  // Finished
  return 0;
}
