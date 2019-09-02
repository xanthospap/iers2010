#ifndef _IERS_1010_HARDISP_
#define _IERS_1010_HARDISP_

/**
 * @file      hardisp.hpp
 * 
 * @version   1.0-1b
 * 
 * @author    xanthos@mail.ntua.gr
 *            danast@mail.ntua.gr
 * 
 * @date      January, 2015
 *
 * @brief     Additional C++ functions and definitions implementing the IERS 
 *            2010 standards; This header file accompanies the hardisp
 *            software.  
 *
 * @note      
 *    -# Original FORTRAN software can be found at 
 *       http://maia.usno.navy.mil/conv2010/software.html
 *    -# The compilation flag <b>QUICK_EXIT</b> can be set to implement a 
 *       quick version but reduced accuracy option. See th readme. file for 
 *       details.
 *    -# For more information, see the iers2010.hpp file at the root directory. 
 * 
 * @attention The FORTRAN subroutines may be updated; see that their C++ 
 *            translations stay updated too.
 * 
 * @todo
 *         -# Translate Fortran's REAL to float and NOT double.
 * 
 * @cite      iers2010
 *
 * @copyright Copyright © 2015 Dionysos Satellite Observatory, <br>
 *            National Technical University of Athens.         <br>
 *            This work is free. You can redistribute it and/or modify it under
 *            the terms of the Do What The Fuck You Want To Public License, 
 *            Version 2, as published by Sam Hocevar. See http://www.wtfpl.net/
 *            for more details.
 * 
 * <b><center><hr>
 * National Technical University of Athens <br>
 *      Dionysos Satellite Observatory     <br>
 *        Higher Geodesy Laboratory        <br>
 *      http://dionysos.survey.ntua.gr
 * </center></b>
 */ 

#include <cmath>
#include <algorithm>

namespace iers2010 {

  namespace hisp {

    // Parameters below set the buffer size for computing the tides recursively 
    
    // nl: the number of harmonics used in the prediction
    constexpr int nl   { 600 };
    
    // nt: this must also be set in the subroutine admint
    constexpr int nt   { 342 };
    
    // ntin: the number of harmonics read in
    constexpr int ntin { 11  };

    int
    toymd(const int*,int*);

    double
    etutc(const double&);

    int
    eval(const double&,const int&,const double*,const double*,
          const double*,double&);

    int
    recurs(double*,const int&,const double*,const int&,const double*,
          double*);

    int
    shells(double*,int*,const int&);

    int
    spline(const int&,const double*,const double*,double*,double*);

    int
    tdfrph(const int*,const int*,double&,double&);
    
    int
    admint(const double*,const int [][6],const double*,double*,
          double*,double*,const int&,int&,const int*);


    inline bool
    leap(int y) noexcept
    {
      int leap = 1 - ((y%4)+3)/4;
      if ( !(y%100) && (y%400) ) leap = 0;
      return leap;
    }

    /**
     * @details This function finds the day number of days before start of 
     *          month m, of year iy, in Gregorian intercalation.
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
    inline int
    mday(int iy, int m) noexcept
    {
      int isleap = leap(iy);
      return (((367*(m-2-12*((m-14)/12)))/12+29) % (365)) + 
        isleap*((9+m)/12);
    }

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
    inline int
    juldat(int y, int m, int d) noexcept
    {
      return (1461*(y+4800+(m-14)/12))/4 + 
          (367*(m-2-12*((m-14)/12)))/12 - 
          (3*((y+4900+(m-14)/12)/100))/4+d-32075;
    }

  } /* end namespace hisp */

}/* end namespace iers2010 */

#endif
