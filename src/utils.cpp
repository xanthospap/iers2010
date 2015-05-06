#include "hardisp.hpp"

/**
 * @details This function converts times expressed in year and
 *          day of year to year-month-day.
 * 
 * @param[in]  it1    Time given in year and day of year (Note 1)
 * @param[out] it2    Time given in year-month-day format (Note 2)
 * @return            Zero
 * 
 * @note
 *     -# The time is split into a year, and the day of the
 *        year,  i.e. it1 is an array of the form [year,  day_of_year].
 *     -# Time is represented as an array of the form [year, month, day]
 *     -# Status:  Class 1 model
 * 
 * Reference: GPS Toolbox:Date/Time conversion algorithms by
 *            Benjamin W. Remondi, see
 *            http://www.ngs.noaa.gov/gps-toolbox/bwr-c.txt
 * 
 */
int iers2010::hisp::toymd (const int* it1, int* it2)
{
  static long month_day[2][13] = {
    {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
    {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}
  };

  long leap  ( (*it1) % 4 == 0 );
  long guess ( it1[1] * 0.032 );
  long more  ( ( (it1[1] - month_day[leap][guess+1]) > 0 ) );
  it2[0]     = it1[0];
  it2[1]     = guess + more + 1;
  it2[2]     = it1[1] - month_day[leap][guess+more];
  
  return 0;
}

/**
 * @details This function determines whether a given integer year is a leap
 *          year.
 * 
 * @param[in]  iy  Year (Note 1)
 * @return         The function returns a boolean variable, denoting if the 
 *                 given year is leap (i.e if the returned value is 0, then the
 *                 year is NOT leap).
 * 
 * @note
 *     -# The year is a Gregorian year (integer).
 *     -# Status: Canonical model.
 * 
 * @version 2009 August 19
 * 
 */
bool iers2010::hisp::leap (const int& iy)
{
  int leap ( 1 - ( (iy % 4) + 3 ) / 4 );

  if ( ( !(iy % 100) ) && (iy % 400) ) 
      return false;

  return leap;
}
