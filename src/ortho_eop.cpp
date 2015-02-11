#include "iers2010.hpp"

#ifdef USE_EXTERNAL_CONSTS
#include "gencon.hpp"
#endif

#ifdef QUICK_EXIT
#include <algorithm>
#endif

/**
 * @details  The purpose of the function is to compute the diurnal and semi-
 *           diurnal variations in Earth Orientation Parameters (x,y, UT1) from
 *           ocean tides.
 *           This function is a translation/wrapper for the fortran ORTHO_EOP
 *           subroutine, found here : 
 *           http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  time  Modified Julian Date
 * @param[out] eop   A vector of size 3, with elements:
 *                   delta_x, in microarcseconds
 *                   delta_y, in microarcseconds
 *                   delta_UT1, in microseconds
 * @return           An integer value which can be:
 *                   Returned Value | Status
 *                   ---------------|----------------------------------------
 *                              -1  | Error; Invalid year
 *                               0  | All ok; a new value has been computed
 *(Only when QUICK_EXIT enabled) 1  | All ok; previous value for angle used.
 *
 * @note 
 *    -# The diurnal and semidiurnal orthoweights fit to the 8 constituents
 *       are listed in Reference 1.
 *    -# Status: Class 1 model
 *
 * @verbatim
 *  Test case:
 *   given input: MJD = 47100D0
 *   expected output: delta_x = -162.8386373279636530D0 microarcseconds
 *                    delta_y = 117.7907525842668974D0 microarcseconds
 *                    delta_UT1 = -23.39092370609808214D0 microseconds
 * @endverbatim
 *
 * @version 2010 March    19
 *
 * @cite iers2010,
 * Ray, R. D., Steinberg, D. J., Chao, B. F., and Cartwright, D. E.,
 * "Diurnal and Semidiurnal Variations in the Earth's Rotation
 *  Rate Induced by Ocean Tides", 1994, Science, 264, pp. 830-832
 *
 */
int iers2010::ortho_eop (const double& time,double* eop)
{
  // check for quick exit
  #ifdef QUICK_EXIT
    static double previous_time = .0e0;
    static double previous_eop[3];
    if ( fabs(previous_time-time)<DATE_MAX_DIFF ) {
      std::copy (previous_eop,previous_eop+3,eop);
      return 1;
    }
  #endif

  static const double orthow[3][12]= {
    {-6.77832e0,-14.86323e0, 0.47884e0,-1.45303e0, 0.16406e0,  0.42030e0,
      0.09398e0, 25.73054e0,-4.77974e0, 0.28080e0, 1.94539e0, -0.73089e0},
    {14.86283e0, -6.77846e0, 1.45234e0, 0.47888e0,-0.42056e0,  0.16469e0,
     15.30276e0, -4.30615e0, 0.07564e0, 2.28321e0,-0.45717e0, -1.62010e0},
    {-1.76335e0,  1.03364e0,-0.27553e0, 0.34569e0,-0.12343e0, -0.10146e0,
     -0.47119e0,  1.28997e0,-0.19336e0, 0.02724e0, 0.08955e0,  0.04726e0}
  };

  // Compute the partials of the tidal variations to the orthoweights
  double h[12];
  iers2010::cnmtx (time,h);
 
  // Compute eop changes
  for (int i=0;i<3;i++) {
    eop[i] = 0e0;
    for (int j=0;j<12;j++) {
      eop[i] += h[j] * orthow[i][j];
    }
  }

  // update quick exit
  #ifdef QUICK_EXIT
    previous_time = time;
    std::copy (eop,eop+3,previous_eop);
  #endif

  // finished
  return 0;
}
