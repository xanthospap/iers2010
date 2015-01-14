#ifndef _IERS_1010
#define _IERS_1010

/**
 * @file      iers2010.hpp
 * 
 * @version   1.0-1b
 * 
 * @author    xanthos@mail.ntua.gr
 *            danast@mail.ntua.gr
 * 
 * @date      December, 2014
 *
 * @brief     C++ functions and definitions implementing the IERS 2010 standards.
 * 
 * @details   
 *
 * @note      
 *         -# Original FORTRAN software can be found at http://maia.usno.navy.mil/conv2010/software.html
 *         -# The compilation flag <b>QUICK_EXIT</b> can be set to implement a quick version but
 *            reduced accuracy option. See th readme. file for details.
 *         -# <b> CLASIFICATION OF FUNCTIONS </b> <br>
 *            In general, Class 1, 2, and 3 models represent physical effects that
 *            act on geodetic parameters while canonical models provide lower-level
 *            representations or basic computations that are used by Class 1, 2, or
 *            3 models.
 *            Class 1 models are those recommended to be used a priori in the
 *            reduction of raw space geodetic data in order to determine
 *            geodetic parameter estimates.
 *            Class 2 models are those that eliminate an observational
 *            singularity and are purely conventional in nature.
 *            Class 3 models are those that are not required as either Class 1 or 2.
 *            Canonical models are accepted as is and cannot be classified as
 *            a Class 1, 2, or 3 model.
 * 
 * @pre       This library can be used as stand-alone. However, a user can collect all
 *            used constants in a seperate file called gencon.hpp and use this file to
 *            declare / define the constants. In this case, the compilation flag
 *            <b>USE_EXTERNAL_CONSTS</b> should be used.
 * 
 * @attention The FORTRAN subroutines may ba updated; see that their C++ translations stay
 *            updated too.
 *
 * @todo
 *        -#  How can we best implement the FORTRAN MOD() fuction in C++ ? Some subroutines, use
 *            both MOD() and DMOD(). We need a translation rule! Affects:
 *            Function    | Comment
 *            ------------|---------------------------------------------------
 *            pmsdnut2    | has both MOD and MOD
 *            utlibr      | MOD
 *            fundarg     | MOD
 *            step2diu    | MOD
 *            arg2        | DMOD
 * 
 *        -#  Provide a test / driver program with differences from the official FORTRAN release,
 *            using the test cases in the modules.
 *        -#  Input time parameters vary (TT, UTC,  Julian dates,  calendar dates,  etc...),. Need
 *           to specify a unique pattern for this.
 * 
 * @cite      iers2010
 *
 * @copyright No COPYRIGHT whatsoever.
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
  
  #ifdef QUICK_EXIT
    #define DATE_MAX_DIFF 1.0e-06
  #endif
  
  /** @brief Function to compute the diurnal lunisolar effect on polar motion.
   */
  int pmsdnut2 (const double&,double*);
  
  /** @brief Function to compute the subdiurnal librations in UT1.
   */
  int utlibr (const double&,double&,double&);
  
  /** @brief  Function to compute the lunisolar fundamental arguments from the model by Simon et al. (1994).
   */
  int fundarg (const double&,double&,double&,double&,double&,double&);
  
  /** @brief  Computes corrections to the coordinates of the CIP to account for Free Core Nutation.
   */
  int fcnnut (const double&,double&,double&,double&,double&);
  
  /** @brief  Computes the angular argument which depends on time for 11 tidal argument calculations.
  */
  int arg2 (const int&,const double&,double*);
  
  /** @brief Computes tidal corrections of station displacements caused by lunar and solar gravitational attraction.
   */
  int dehanttideinel (const double*,const double*,const double*,const double&,const double&,double*);
  int dehanttideinel (const double*,const double*,const double*,const int&,const int&,const int&,const double&,double*);
  
  /** @brief Computes the effects of zonal Earth tides on the rotation of the Earth.
   */
  int rg_zont2 (const double&,double&,double&,double&);
};

#endif
