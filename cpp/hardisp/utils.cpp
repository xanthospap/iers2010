/**
 * @details This function converts a Gregorian date to a Julian date.
 * 
 * @param[in] it A Gregorian date (Note 1)
 * @return       The corresponding Julian date (Note 2)
 * 
 * @note
 *     -# The format of the Gregorian date should be yyyy-mm-dd (i.e a
 *       3-dimensional array of ints)
 *     -# The date is valid for all positive values.
 *     -# Status: Canonical model
 * 
 * @version 2009 August 19
 * 
 * @cite iers2010, 
 * Explanatory Supplement American Ephemeris & Nautical Almanac
 * (cf Comm CACM, 11, 657 (1968) and 15, 918 (1972)), p. 604
 *
 */
inline int juldat (const int* it)
{
  return (1461*(it[0]+4800+(it[1]-14)/12))/4
         + (367*(it[1]-2-12*((it[1]-14)/12)))/12
         - (3*((it[0]+4900+(it[1]-14)/12)/100))/4+it[2]-32075;
}

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
 * @cite GPS Toolbox:Date/Time conversion algorithms by Benjamin W. Remondi
 * 
 */
int toymd (const int* it1, int* it2)
{
  static long month_day[2][13] = {
    {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
    {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}
  };

  long leap  = ( (*it1) % 4 == 0 );
  long guess = it1[1] * 0.032;
  long more  = ( ( it1[1] - month_day[leap][guess+1] ) > 0 );
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
 * @return         The function returns a boolean variable, denoting if the given
 *                 year is leap (i.e if the returned value is 0, then the year is NOT
 *                 leap).
 * 
 * 
 * @note
 *     -# The year is a Gregorian year (integer).
 *     -# Status: Canonical model.
 * 
 * @version 2009 August 19
 * 
 */
bool leap (const int& iy)
{
  int leap = 1 - ( (iy % 4) + 3 ) / 4;

  if ( ( !(iy % 100) ) && (iy % 400) ) return false;

  return leap;
}

/**
 * @details This function finds the day number of days before start of month m,
 *          of year iy, in Gregorian intercalation.
 * 
 * @param[in]  iy Given year
 * @param[in]  m  Given month
 * @return        Day number of day before start of a month
 * 
 * @note Status: Canonical model
 * 
 * @version 2009 July  29
 * 
 */
inline int mday (const int& iy, const int& m)
{
  int leap = leap (iy);
  return (((367*(m-2-12*((m-14)/12)))/12+29) %  365) + leap*((9+m)/12;
}