#include "iers2010.hpp"

/**
 * @details  This function computes the global total FCULa mapping function 
 *           (Mendes et al. 2002). It is dependent on latitude, height, and 
 *           surface temperature.
 *           This function is a translation/wrapper for the fortran FCUL_A
 *           subroutine, found here : 
 *           http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  dlat  Latitude given in degrees (North Latitude)
 * @param[in]  dhgt  Height given in meters (mean sea level)
 * @param[in]  t     Surface temperature given in Kelvin
 * @param[in]  elev  Elevation angle in degrees
 * @param[out] fcul  Mapping function to scale total delay (Note 1)
 * @return           An integer, always zero
 * 
 * @note
 *    -# These coefficients are based on a LS adjustment of 87766 (cleaned)
 *    set of traces, based on Ciddor routines to compute refractivity,
 *    according to IUGG recommendations (1999).
 *    -# Status: Class 1 model
 *
 * @verbatim
 *    given input: LATITUDE = 30.67166667D0 degrees (McDonald Observatory)
 *                  HEIGHT_M = 2075D0 meters (mean sea level)
 *                  T_K      = 300.15D0 Kelvin (August 12, 2009)
 *                  ELEV_DEG = 15D0 degrees (See Mendes et al.)
 *     expected output: FCUL_A = 3.800243667312344087D0
 * @endverbatim
 * 
 * @version 2009 August 13
 *
 * @cite iers2010
 * Mendes, V.B., G. Prates, E.C. Pavlis, D.E. Pavlis, 
 * and R.B. Langley (2002). "Improved mapping functions for
 * atmospheric refraction correction in SLR", Geophysical 
 * Res. Lett., 29(10), 1414, doi:10.1029/2001GL014394, 2002
 * 
 */
int iers2010::fcul_a (const double& dlat,const double& dhgt,const double& t,
    const double& elev,double& fcul)
{
    #ifdef USE_EXTERNAL_CONSTS
        constexpr double PI   (DPI);
    #else
        constexpr double PI   (3.14159265358979323846e0);
    #endif
    
    // Convert elevation angle to radians
    double epsilon  ( elev * (PI/180e0) );
    double sine     ( sin (epsilon) );
    // Convert temperature to Celsius
    double t_c      ( t - 273.15e0 );
    double cosphi   ( cos (dlat * (PI/180e0)) );

    // Define coefficients used in the model
    double a10 (  0.121008e-02 );
    double a11 (  0.17295e-05  );
    double a12 (  0.3191e-04   );
    double a13 ( -0.18478e-07  );

    double a20 (  0.304965e-02 );
    double a21 (  0.2346e-05   );
    double a22 ( -0.1035e-03   );
    double a23 ( -0.1856e-07   );
    
    double a30 (  0.68777e-01 );
    double a31 (  0.1972e-04  );
    double a32 ( -0.3458e-02  );
    double a33 (  0.1060e-06  );

    // a, b, and c in Marini continued fraction (Eq. 5)
    double a1 ( a10+a11*t_c+a12*cosphi+a13*dhgt );
    double a2 ( a20+a21*t_c+a22*cosphi+a23*dhgt );
    double a3 ( a30+a31*t_c+a32*cosphi+a33*dhgt );

    // numerator in continued fraction
    double map_zen ( (1.0e0 + a1/(1.0e0 + a2/(1.0e0+a3))) );

    // result
    fcul = map_zen/(sine+a1/(sine+a2/(sine+a3)));

    // return
    return 0;
}
