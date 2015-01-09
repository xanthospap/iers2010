#include "hardisp.hpp"

/**
 * @details The purpose of the function is to compute the difference, delta,
 *          between Epheremis Time (ET) and Coordinated Universal Time (UTC).
 * 
 * @param[in]  year  Decimal year (Note 1)
 * @return           Delta ET - UTC (Note 2)
 * 
 * @note
 *     -# This function is valid only from 1700.-until next leap second.
 *        Currently, this is up to 2012.5.
 *     -# The expression used in given in seconds.
 *     -# Leap second table in GAMIT UTC (and UT) is the time most 
 *        often used (e.g. in time signals)
 *     -# Status: Canonical model
 * 
 * @warning \b IMPORTANT <br>
 *     A new version of this routine must be
 *     produced whenever a new leap second is
 *     announced.  There are three items to
 *     change on each such occasion:<br>
 *     1) Update the nstep variable<br>
 *     2) Update the arrays st and si<br>
 *     3) Change date of latest leap second<br>
 *     <b>Latest leap second:  2012 June 30</b>
 * 
 * Test case:
 *     given input: year = 2007.0 
 *
 *     expected output: delta = 65.1840000000000 seconds
 *
 *     given input: year = 2013.0 
 *
 *     expected output: delta = 67.1840000000000 seconds
 *
 * @version 2012 March 13
 * 
 */
  double iers2010::hisp::etutc (const double& year)
  {
    const int nstep = 25;
    
    // si gives amount of step, at the times given in st
    const double si = 25e0*1e0;
    
    static const double st[] = {
      1972.5,1973.0,1974.0,1975.0,1976.0,1977.0,1978.0,
      1979.0,1980.0,1981.5,1982.5,1983.5,1985.5,1988.0,
      1990.0,1991.0,1992.5,1993.5,1994.5,1996.0,1997.5,
      1999.0,2006.0,2009.0,2012.5
    };
    
    static const double d[142] = {
       5.15e0, 4.64e0, 5.36e0, 3.49e0, 3.27e0, 2.45e0, 4.03e0, 1.76e0, 3.30e0,
       1.00e0, 2.42e0, 0.94e0, 2.31e0, 2.27e0,-0.22e0, 0.03e0,-0.05e0,-0.06e0,-0.57e0,
       0.03e0,-0.47e0, 0.98e0,-0.86e0, 2.45e0, 0.22e0, 0.37e0, 2.79e0, 1.20e0, 3.52e0,
       1.17e0, 2.67e0, 3.06e0, 2.66e0, 2.97e0, 3.28e0, 3.31e0, 3.33e0, 3.23e0, 3.60e0,
       3.52e0, 4.27e0, 2.68e0, 2.75e0, 2.67e0, 1.94e0, 1.39e0, 1.66e0, 0.88e0, 0.33e0,
      -0.17e0,-1.88e0,-3.43e0,-4.05e0,-5.77e0,-7.06e0,-7.36e0,-7.67e0,-7.64e0,-7.93e0,
      -7.82e0,-8.35e0,-7.91e0,-8.03e0,-9.14e0,-8.18e0,-7.88e0,-7.62e0,-7.17e0,-8.14e0,
      -7.59e0,-7.17e0,-7.94e0,-8.23e0,-7.88e0,-7.68e0,-6.94e0,-6.89e0,-7.11e0,-5.87e0,
      -5.04e0,-3.90e0,-2.87e0,-0.58e0, 0.71e0, 1.80e0, 3.08e0, 4.63e0, 5.86e0, 7.21e0,
       8.58e0,10.50e0,12.10e0,12.49e0,14.41e0,15.59e0,15.81e0,17.52e0,19.01e0,18.39e0,
      19.55e0,20.36e0,21.01e0,21.81e0,21.76e0,22.35e0,22.68e0,22.94e0,22.93e0,22.69e0,
      22.94e0,23.20e0,23.31e0,23.63e0,23.47e0,23.68e0,23.62e0,23.53e0,23.59e0,23.99e0,
      23.80e0,24.20e0,24.99e0,24.97e0,25.72e0,26.21e0,26.37e0,26.89e0,27.68e0,28.13e0,
      28.94e0,29.42e0,29.66e0,30.29e0,30.96e0,31.09e0,31.59e0,32.06e0,31.82e0,32.69e0,
      33.05e0,33.16e0,33.59
    };
    
    static const double tx[39] = {
      61.5e0,
      62.e0     ,62.5e0     ,63.e0      ,63.5e0     ,64.e0      ,64.5e0     ,65.e0   ,
      65.5e0    ,66.e0      ,66.5e0     ,67.e0      ,67.5e0     ,68.e0      ,68.25e0 ,
      68.5e0    ,68.75e0    ,69.e0      ,69.25e0    ,69.5e0     ,69.75e0    ,70.e0   ,
      70.25e0   ,70.5e0     ,70.75e0    ,71.e0      ,71.085e0   ,71.162e0   ,71.247e0,
      71.329e0  ,71.414e0   ,71.496e0   ,71.581e0   ,71.666e0   ,71.748e0   ,71.833e0,
      71.915e0  ,71.999e0   ,72.0e0
    };
    
    static const double ty[39] = {
      33.59e0,
      34.032e0  ,34.235e0   ,34.441e0   ,34.644e0   ,34.95e0    ,35.286e0   ,35.725e0,
      36.16e0   ,36.498e0   ,36.968e0   ,37.444e0   ,37.913e0   ,38.39e0    ,38.526e0,
      38.76e0   ,39.e0      ,39.238e0   ,39.472e0   ,39.707e0   ,39.946e0   ,40.185e0,
      40.42e0   ,40.654e0   ,40.892e0   ,41.131e0   ,41.211e0   ,41.284e0   ,41.364e0,
      41.442e0  ,41.522e0   ,41.600e0   ,41.680e0   ,41.761e0   ,41.838e0   ,41.919e0,
      41.996e0  ,42.184e0   ,42.184e0
    };
    
    // For the oldest epochs, use approximations
    if ( year < 1700 ) 
      return /*delta =*/ .0e0;
      
    else if ( year < 1785.0)
      return /*delta =*/ (year-1750.0e0)/5.0e0;
      
    else if ( year < 1820.5 )
      return /*delta =*/ 6.0e0;
      
    // For 1820.5 to 1961.5, data is spaced at yearly intervals
    else if ( year < 1961.5 ) {
      double n = year - 1819.5e0;
      double frac  = year - (1819.5e0 + n);
      double delta = (d[n] - d[n-1])*frac + d[n-1];
      return delta;
    }

    // For 1961.5 to 1972.0, interpolate between unequispaced data
    else if ( year < 1972.0 ) {
      for (int i = 0;i<38;i++) {
        if ( (year-1900.0) == tx[i] ) {
          return /*delta =*/ty[i];
        }
        if ( (year-1900.0) < tx[i] ) {
          delta = ty[i-1] + (ty[i]-ty[i-1])*
                  ((year-1900.0e0-tx[i-1])/(tx[i]-tx[i-1]));
          return delta;
        }
      }
    }

  /*--------------------------------------------------------------------------*
   *   after 1972 et-utc has only step offsets. st is the array of step times,*
   *   and si is the step sizes (an added second is +1.)                      *
   *--------------------------------------------------------------------------*/
  double delta = 42.184e0;
  for (int i=0;i<nstep;i++) {
    if (year >= st[i]) delta += si[i];
    if (year < st[i]) return delta;
  }
  
  return delta;
}