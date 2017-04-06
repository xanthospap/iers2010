#include "iers2010.hpp"

/**
 * @details Gregorian Calendar to Julian Date. This function is closely related
 *          to the SOFA function iauCal2jd. 
 *
 * @param[in]     iy   Year
 * @param[in]     im   Month
 * @param[in]     id   Day
 * @param[out]    djm0 MJD zero-point: always 2400000.5
 * @param[out]    djm  Modified Julian Date for 0 hrs
 * @return             An \c integer denoting the exit status
 * Returned Int  | Status
 * ------------- | -------------
 *  0            | OK
 * -1            | Bad year (Note 3: JD not computed)
 * -2            | Bad month (JD not computed)
 * -3            | Bad day (JD computed)
 *
 * @note
 *  -# The algorithm used is valid from -4800 March 1, but this
 *     implementation rejects dates before -4799 January 1.
 *  -# The Julian Date is returned in two pieces, in the usual SOFA
 *     manner, which is designed to preserve time resolution.  The
 *     Julian Date is available as a single number by adding djm0 and
 *     djm.
 *  -# In early eras the conversion is from the "Proleptic Gregorian
 *     Calendar";  no account is taken of the date(s) of adoption of
 *     the Gregorian Calendar, nor is the AD/BC numbering convention
 *     observed.
 *
 * @cite esaa , Section 12.92 (p604).
 *
 * @version 07.11.2014 (release  2016-05-03)
 */
int
iers2010::dhtide::cal2jd(int iy, int im, int id, double& djm0,double& djm)
{

    // Earliest year allowed (4800BC)
    constexpr int IYMIN { -4799 };

    // Month lengths in days
    /*static*/ constexpr int mtab[] = {
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
    };

    // Preset status.
    int j { 0 };

    // Validate year and month
    if (iy < IYMIN) { 
        return -1;
    }
    if (im < 1 || im > 12) {
        return -2;
    }

    // If February in a leap year, 1, otherwise 0.
    int ly { ( (im == 2) && !(iy%4) && (iy%100 || !(iy%400)) ) };

    // Validate day, taking into account leap years.
    if ( (id < 1) || (id > (mtab[im-1] + ly)) ) {
        j = -3;
    }

    // Return result.
    int  my    { (im - 14) / 12 };
    long iypmy { (long)(iy + my) };
    djm0 = 2400000.5e0; // or DJM0
    djm  = (double)((1461L * (iypmy + 4800L)) / 4L
                  + (367L  * (long)(im - 2 - 12*my)) / 12L
                  - (3L    * ((iypmy + 4900L) / 100L)) / 4L
                  + (long)id - 2432076L);
    
    // Return status
    return j;
}
