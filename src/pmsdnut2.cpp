#include "iers2010.hpp"

#ifdef USE_EXTERNAL_CONSTS
  #include "gencon.hpp"
#endif

/**
 * @details This function evaluates the model of polar motion for
 *          a nonrigid Earth due to tidal gravitation. This polar motion
 *          is equivalent to the so-called "subdiurnal nutation." The model
 *          is a sum of a first order polynomial and 25 trigonometric terms
 *          (15 long periodic and 10 quasi diurnal) with coefficients given
 *          in Table 5.1a of the IERS Conventions (2010).
 *          This function is a translation/wrapper for the fortran PMSDNUT2
 *          subroutine, found here : 
 *          http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  rmjd Time expressed as modified Julian date
 * @param[out] dx   The x component of polar motion expressed in 
 *                  microarcseconds.
 * @param[out] dy   The y component of polar motion expressed in 
 *                  microarcseconds.
 * @return          An integer value, always 0.
 *
 * @note    Status:  Class 1 model
 * 
 * @warning In the present version this subroutine neglects the linear trend
 *          and the long periodic terms of the expansion, for the reasons 
 *          explained in Section 5.5.1.1 of the IERS Conventions (2010). If 
 *          the full expansion is needed, set the parameter iband to 0 instead
 *          of 1, that is, replace the statement
 *          PARAMETER ( iband = 1 )
 *          to  PARAMETER ( iband = 0 )
 * 
 * @version 13.10.2011
 * 
 * @cite iers2010
 * 
 */
int
iers2010::pmsdnut2(double rmjd, double& dx, double& dy) 
{
    /*
     *         ----------------------------
     *           D E F I N I T I O N S
     *         ----------------------------
     *
     *  iband  - parameter defining the range of periods for the terms which
     *           are included in computations; if equal to 1 only the quasi 
     *           diurnal terms are computed, otherwise the full model
     *  iarg   - array defining for each of the 25 trigonometric terms a set
     *           of 6 integer multipliers of the fundamental angular arguments
     *  arg    - vector of the following 6 fundamental arguments used to
     *           compute the angular argument of the trigonometric functions
     *           arg(1:6) = [ GMST+pi, el, elp, f, d, om ]; this vector is
     *           evaluated by the subroutine FUNDARG which is called as an 
     *           external subroutine.  Originally evaluated by the subroutine
     *           PMARGS. 
     *  period - array of periods of the trigonometric terms of expansion, in
     *           mean solar days; only for a check - not used in computations
     *  xs, xc - sine and cosine coefficients of the x coordinate of the pole,
     *           in microarcseconds
     *  ys, yc - sine and cosine coefficients of the y coordinate of the pole,
     *           in microarcseconds
     *  angle  - angular argument of the trigonometric functions
     *           angle = Sum(i=1:6) iarg(i,j)*arg(i), for j=1,25
     * 
     */
    int    iband = 1;
    double fargs[6]; /* fundamental arguments: GMST+pi, l, lp, f, d, om */

    // Set constants
    #ifdef USE_EXTERNAL_CONSTS
        constexpr double RMJD0   (DJM00);     // Modified Julian date of J2000
        constexpr double PI      (DPI);
        constexpr double TWOPI   (D2PI); 
        constexpr double RAD2SEC (DRAD2SEC);  // Radians to seconds
    #else
        // Modified Julian date of J2000
        constexpr double RMJD0   (51544.5e0);
        constexpr double PI      (3.141592653589793238462643e0);
        constexpr double TWOPI   (6.283185307179586476925287e0);
        // Radians to seconds
        constexpr double RAD2SEC (86400e0 / TWOPI);
    #endif

    // Coefficients of the long and quasi diurnal periodic terms in polar motion
    constexpr struct {
        int    iarg[6];
        double per, xs, xc, ys, yc;
    } x[] = {
        //  Coefficients of the long periodic terms in polar motion
        //+ Source: IERS Conventions (2010), Table 5.1a
        { { 0,  0, 0,  0,  0, -1 }, 6798.3837,    0.0,   0.6,   -0.1,   -0.1 },
        { { 0, -1, 0,  1,  0,  2 }, 6159.1355,    1.5,   0.0,   -0.2,    0.1 },
        { { 0, -1, 0,  1,  0,  1 }, 3231.4956,  -28.5,  -0.2,    3.4,   -3.9 },
        { { 0, -1, 0,  1,  0,  0 }, 2190.3501,   -4.7,  -0.1,    0.6,   -0.9 },
        { { 0,  1, 1, -1,  0,  0 }, 438.35990,   -0.7,   0.2,   -0.2,   -0.7 },
        { { 0,  1, 1, -1,  0, -1 }, 411.80661,    1.0,   0.3,   -0.3,    1.0 },
        { { 0,  0, 0,  1, -1,  1 }, 365.24219,    1.2,   0.2,   -0.2,    1.4 },
        { { 0,  1, 0,  1, -2,  1 }, 193.55971,    1.3,   0.4,   -0.2,    2.9 },
        { { 0,  0, 0,  1,  0,  2 }, 27.431826,   -0.1,  -0.2,    0.0,   -1.7 },
        { { 0,  0, 0,  1,  0,  1 }, 27.321582,    0.9,   4.0,   -0.1,   32.4 },
        { { 0,  0, 0,  1,  0,  0 }, 27.212221,    0.1,   0.6,    0.0,    5.1 },
        { { 0, -1, 0,  1,  2,  1 }, 14.698136,    0.0,   0.1,    0.0,    0.6 },
        { { 0,  1, 0,  1,  0,  1 }, 13.718786,   -0.1,   0.3,    0.0,    2.7 },
        { { 0,  0, 0,  3,  0,  3 }, 9.1071941,   -0.1,   0.1,    0.0,    0.9 },
        { { 0,  0, 0,  3,  0,  2 }, 9.0950103,   -0.1,   0.1,    0.0,    0.6 },
        //  Coefficients of the quasi diurnal terms in polar motion
        //+ Source: IERS Conventions (2010), Table 5.1a
        { { 1, -1, 0, -2,  0, -1 }, 1.1196992,   -0.4,   0.3,   -0.3,   -0.4 },
        { { 1, -1, 0, -2,  0, -2 }, 1.1195149,   -2.3,   1.3,   -1.3,   -2.3 },
        { { 1,  1, 0, -2, -2, -2 }, 1.1134606,   -0.4,   0.3,   -0.3,   -0.4 },
        { { 1,  0, 0, -2,  0, -1 }, 1.0759762,   -2.1,   1.2,   -1.2,   -2.1 },
        { { 1,  0, 0, -2,  0, -2 }, 1.0758059,  -11.4,   6.5,   -6.5,  -11.4 },
        { { 1, -1, 0,  0,  0,  0 }, 1.0347187,    0.8,  -0.5,    0.5,    0.8 },
        { { 1,  0, 0, -2,  2, -2 }, 1.0027454,   -4.8,   2.7,   -2.7,   -4.8 },
        { { 1,  0, 0,  0,  0,  0 }, 0.9972696,   14.3,  -8.2,    8.2,   14.3 },
        { { 1,  0, 0,  0,  0, -1 }, 0.9971233,    1.9,  -1.1,    1.1,    1.9 },
        { { 1,  1, 0,  0,  0,  0 }, 0.9624365,    0.8,  -0.4,    0.4,    0.8 }
    };
    constexpr int M { sizeof(x) / sizeof(x[0]) };
    static_assert( M == 25, "Invalid quasi diurnal terms in pmsdnut2." );

    //  Rate of secular polar motion, in microarcseconds per year
    //+ Source: IERS Conventions (2010), Table 5.1a
    constexpr double xrate { -3.8e0 },
                     yrate { -4.3e0 };

    //  Compute the periodical part of the model
    //+ Coordinates of the pole are set to zero first
    dx = dy = .0e0;

    //  Evaluate the vector of the fundamental arguments
    //+ arg(0:5) = [ GMST+pi, el, elp, f, d, om ]
    //+ at t = rmjd

    // Convert the input epoch to Julian centuries of TDB since J2000
    double t { (rmjd-RMJD0) / 36525e0 };

    // Compute GMST + pi
    double gmst { std::fmod(67310.54841e0 +
                    t*((8640184.812866e0 + 3155760000e0) +
                    t*(0.093104e0 +
                    t*(-0.0000062 ))), 86400e0) };

    // Fundamental arguments
    iers2010::fundarg(t, fargs+1);
    fargs[0] = gmst / RAD2SEC + PI;
    fargs[0] = std::fmod(fargs[0], TWOPI);

    int jstart { (iband == 1) ? 15 : 0 };

    double angle, sina, cosa;
    for (int j=jstart; j<M; j++) {
        //  For the j-th term of the trigonometric expansion, compute the angular
        //+ argument angle of sine and cosine functions as a linear integer
        //+ combination of the 6 fundamental arguments
        angle = .0e0;
        for (int i=0; i<6; i++) {
            angle += ( double(x[j].iarg[i]) * fargs[i] );
        }
        angle = std::fmod(angle, TWOPI);
        // Compute contribution from the j-th term to the polar motion coordinates
        sina = std::sin(angle);
        cosa = std::cos(angle);
        dx += x[j].xs * sina + x[j].xc * cosa;
        dy += x[j].ys * sina + x[j].yc * cosa;
    }

    if (iband != 1) {
        // Add the secular term of the model
        dx += xrate * (rmjd-RMJD0) / 365.25e0;
        dy += yrate * (rmjd-RMJD0) / 365.25e0;
    }

    //  Finished
    return 0;
}
