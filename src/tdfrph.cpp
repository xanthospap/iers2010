#include "iers2010.hpp"

/**
 * \details This function converts times expressed in year and
 *          day of year to year-month-day.
 * 
 * \param[in]  iy     The year.
 * \param[in]  doy    The day of year.
 * \param[out] month  The month.
 * \param[out] dom    The day of month.
 * \return            Zero
 * 
 * \note
 *     -# Status:  Class 1 model
 * 
 * Reference: GPS Toolbox:Date/Time conversion algorithms by
 *            Benjamin W. Remondi, see
 *            http://www.ngs.noaa.gov/gps-toolbox/bwr-c.txt
 * 
 */
int
toymd(int iy, int doy, int& month, int& dom)
{
    static long month_day[2][13] = {
        {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
        {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}
    };

    long leap  { iy % 4 == 0 };
    long guess ( (long)doy * 0.032 );
    long more  { ( (doy - month_day[leap][guess+1]) > 0 ) };
    month      = static_cast<int>(guess + more + 1);
    dom        = static_cast<int>(doy - month_day[leap][guess+more]);

    return 0;
}

/**
 * \details This function converts a Gregorian date to a Julian date.
 * 
 * \param[in]  iy     The year.
 * \param[in]  month  The month.
 * \param[in]  dom    The day of month.
 * \return            The corresponding Julian date (as long)
 */
inline long
juldat(int iy, int im, int id)
{
    return (1461L*(iy+4800L+(im-14L)/12L)) / 4L
          +(367L*(im-2L-12L*((im-14L)/12L))) / 12L
          -(3L*((iy+4900L+(im-14L)/12L)/100L)) / 4L + id - 32075L;
}

/**
 * \details This function determines whether a given integer year is a leap year.
 *
 * \param[in]   iy     The (input) year
 * \return      True if year is leap, false otherwise.
 */
inline bool
leap(int iy)
{
    int leap { 1-((iy%4)+3) / 4 };
    if ( !(iy % 100) && (iy % 400) ) { leap = 0; }
    return (bool)leap;
}

/**
 * \details This subroutine returns the frequency and phase of a tidal
 *          constituent when its Doodson number is given as input. 
 * 
 * \param[in]  idood Doodson number of a tidal constituent (6-element integer array)
 * \param[in]  itm   Date as integer array in UTC. The format should be:
 *                   [year, day_of_year, hours, minutes, seconds].
 * \param[out] freq  Frequency of a tidal constituent
 * \param[out] phase Phase of a tidal constituent (Note 1)
 * \return           Always 0.
 * 
 * \note
 *     -# The phases must be decreased by 90 degrees if the sum of the order 
 *        and the species number is odd (as for the 2nd degree diurnals, and 
 *        3rd degree low frequency and semidiurnals).
 *        These phases may need further adjustment to allow for the spherical
 *        harmonic normalization used; e.g. for that used for the potential
 *        by Cartwright and Tayler, 180 degrees must be added for (species,
 *        order) = (1,2), (1,3), or (3,3).
 *     -# Status:  Class 1 model
 * 
 * \version 11.09.2013
 *
 * \cite iers2010
 * 
 * @todo There is a bug here! The line: 
 * 'double delta  ( iers2010::hisp::etutc (year) );'
 * returns a delta != 0. The FORTRAN routine however gives delta=0.
 * WTF ??
 */

int
iers2010::hisp::tdfrph(const int idood[6], const int itm[5], double& freq,
    double& phase)
{
    static int last_date[] = {0, 0, 0, 0, 0};
    static double  d[6];
    static double dd[6];
    
    /*------------------------------------------------------------------------
     *  Test to see if time has changed; if so, set the phases and frequencies
     *  for each of the Doodson arguments
     *------------------------------------------------------------------------*/
    // bool reset_date { false };
    if (   itm[0] != last_date[0]
        || itm[1] != last_date[1] 
        || itm[2] != last_date[2] 
        || itm[3] != last_date[3] 
        || itm[4] != last_date[4] )
     {
        // reset_date = true;

        // Convert times to Julian days (UT) then to Julian centuries 
        // from J2000.0 (ET)
        int month, day_of_month;
        toymd(itm[0], itm[1], month, day_of_month);
        int    jd     { (int)juldat(itm[0], month, day_of_month) };
        double dayfr  { itm[2]/24.0e0 + itm[3]/(24.0*60.0e0) + itm[4]/84600.0e0 };
        double year   { itm[0]+(itm[1]+dayfr)/( 365.0e0 + (double)leap(itm[0])) };
        double delta  { iers2010::hisp::etutc(year) };
        double djd    { (double)jd - 0.5e0 + dayfr  };
        // delta = .0e0;
        double t      { (djd - 2451545.0e0 + delta/86400.0e0)/36525.0e0 };

         // IERS expressions for the Delaunay arguments, in degrees
        double f1 {     134.9634025100e0 +
                 t*( 477198.8675605000e0 +
                 t*(      0.0088553333e0 +
                 t*(      0.0000143431e0 +
                 t*(     -0.0000000680e0 ))))
                 };

        double f2 {   357.5291091806e0 +
                t*( 35999.0502911389e0 +
                t*(    -0.0001536667e0 +
                t*(     0.0000000378e0 +
                t*(    -0.0000000032e0 ))))
                };

        double f3 {     93.2720906200e0 + 
                t*( 483202.0174577222e0 +
                t*(     -0.0035420000e0 +
                t*(     -0.0000002881e0 +
                t*(      0.0000000012e0 ))))
                };

        double f4 {    297.8501954694e0 +
                t*( 445267.1114469445e0 +
                t*(     -0.0017696111e0 +
                t*(      0.0000018314e0 +
                t*(     -0.0000000088e0 ))))
                };

        double f5 {   125.0445550100e0 +
                t*( -1934.1362619722e0 +
                t*(     0.0020756111e0 +
                t*(     0.0000021394e0 +
                t*(    -0.0000000165e0 ))))
                };
            
        // Convert to Doodson (Darwin) variables
        d[0] = 360.0e0*dayfr - f4;
        d[1] = f3 + f5;
        d[2] = d[1] - f4;
        d[3] = d[1] - f1;
        d[4] = -f5;
        d[5] = d[2] - f2;
        
        //  Find frequencies of Delauney variables (in cycles/day), and from
        //+ these the same for the Doodson arguments
        double fd1 {  0.0362916471e0 + 0.0000000013e0*t };
        double fd2 {  0.0027377786e0 };
        double fd3 {  0.0367481951e0 - 0.0000000005e0*t };
        double fd4 {  0.0338631920e0 - 0.0000000003e0*t };
        double fd5 { -0.0001470938e0 + 0.0000000003e0*t };
        dd[0] = 1.0e0 - fd4;
        dd[1] = fd3 + fd5;
        dd[2] = dd[1] - fd4;
        dd[3] = dd[1] - fd1;
        dd[4] = -fd5;
        dd[5] = dd[2] - fd2;
    
        // copy the just used date to the static one for next use of function
        last_date[0] = itm[0];
        last_date[1] = itm[1];
        last_date[2] = itm[2];
        last_date[3] = itm[3];
        last_date[4] = itm[4];

    } // End of intialization (likely to be called only once)
    
    //  Compute phase and frequency of the given tidal constituent
    freq  = 0.0e0;
    phase = 0.0e0;
    for (int i=0; i<6; i++) {
        freq  += ((double)idood[i]) * dd[i];
        phase += ((double)idood[i]) * d[i];
    }
    
    // Adjust phases so that they fall in the positive range 0 to 360
    phase = std::fmod(phase, 360.0e0);
    if ( phase < 0.0e0 ) { phase += 360.0e0; }
    
    // Finished
    return 0;
}
