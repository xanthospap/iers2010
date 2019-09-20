#ifndef _IERS_1010_HARDISP_
#define _IERS_1010_HARDISP_
#include "ggdatetime/dtcalendar.hpp"

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

namespace iers2010 {

  namespace hisp {

    // Parameters below set the buffer size for computing the tides recursively 
    
    // nl: the number of harmonics used in the prediction
    constexpr int nl   { 600 };
    
    // nt: this must also be set in the subroutine admint
    constexpr int nt   { 342 };
    
    // ntin: the number of harmonics read in
    constexpr int ntin { 11  };
    
    double
    eval(double, int, const double*, const double*, const double*);

    int
    recurs(double*, int, const double*, int, const double*, double*);

    int
    shells(double*, int*, int) noexcept;

    int
    spline(int, const double*, const double*, double*, double*);

    int
    tdfrph(const int*, ngpt::datetime<ngpt::seconds>, double&, double&);
    
    int
    admint(const double*, const double*, ngpt::datetime<ngpt::seconds>, double*,
      double*, double*, int, int&);

    int
    read_hardisp_args(double tamp[3][ntin], double tph[3][ntin],
      const char* filename=nullptr);

    int
    hardisp_impl(int, double, double tamp[3][ntin], double tph[3][ntin], 
      ngpt::datetime<ngpt::seconds> epoch);
  } /* end namespace hisp */

}/* end namespace iers2010 */

#endif
