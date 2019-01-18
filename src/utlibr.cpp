#include "iers2010.hpp"

#ifdef USE_EXTERNAL_CONSTS
  #include "gencon.hpp"
#endif

/**
 * @details This function evaluates the model of subdiurnal libration
 *          in the axial component of rotation, expressed by UT1 and LOD.
 *          This effect is due to the influence of tidal gravitation on the
 *          departures of the Earth's mass distribution from the rotational
 *          symmetry, expressed by the non-zonal components of geopotential.
 *          The amplitudes have been computed for an elastic Earth with liquid
 *          core. The adopted truncation level is 0.033 microseconds in UT1
 *          corresponding to the angular displacement of 0.5 microarcseconds
 *          or to 0.015 mm at the planet surface. With this truncation level
 *          the model contains 11 semidiurnal terms. The coefficients of
 *          the model are given in Table 5.1b of the IERS Conventions (2010).
 *          This function is a translation/wrapper for the fortran UTLIBR
 *          subroutine, found here : 
 *          http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  rmjd Time expressed as modified Julian date
 * @param[out] dut1 Incremental UT1 in microseconds
 * @param[out] dlod Incremental LOD in microseconds per day
 * @return          An integer value, always 0. 
 *
 * @note
 *      - Status:  Class 3 model
 * 
 * @version 23.06.2010
 * 
 * @cite iers2010
 * 
 */
int
iers2010::utlibr(double rmjd, double& dut1, double& dlod) 
{
    /*
     *         ----------------------------
     *           D E F I N I T I O N S
     *         ----------------------------
     *  iarg   - array defining for each of the 11 trigonometric terms a set
     *           of 6 integer multipliers of the fundamental angular arguments
     *  arg    - vector of the following 6 fundamental arguments used to
     *           compute the angular argument of the trigonometric functions
     *           arg(1:6) = [ GMST+pi, el, elp, f, d, om ]; this vector is
     *           evaluated by the subroutine FUNDARG which is called as an 
     *           external subroutine.  Originally evaluated by the subroutine
     *           PMARGS. 
     *  period - array of periods of the trigonometric terms of expansion, in
     *           mean solar days; only for a check - not used in computations
     *  dUT1s, dUT1c - sine and cosine coefficients of dUT1, in microseconds
     *  dLODs, dLODc - sine and cosine coefficients of dLOD, in microseconds
     *                 per day
     *  angle  - angular argument of the trigonometric functions
     *           angle = Sum(i=1:6) iarg(i,j)*arg(i), for j=1,11
     */

    // Set constants
    #ifdef USE_EXTERNAL_CONSTS
        // Modified Julian date of J2000
        constexpr double RMJD0   (DJM00);
        constexpr double PI      (DPI);
        constexpr double TWOPI   (D2PI);
        // Radians to seconds
        constexpr double RAD2SEC (DRAD2SEC);
    #else
        // Modified Julian date of J2000
        constexpr double RMJD0   (51544.5e0);
        constexpr double PI      (3.141592653589793238462643e0);
        constexpr double TWOPI   (6.283185307179586476925287e0);
        // Radians to seconds
        constexpr double RAD2SEC (86400e0 / TWOPI);
    #endif

    //  Coefficients of the quasi semidiurnal terms in dUT1, dLOD 
    //+ Source: IERS Conventions (2010), Table 5.1b
    constexpr struct {
        int iarg[6];
        double per, dut1s, dut1c, dlods, dlodc;
    } x[] = {
        { { 2, -2,  0, -2,  0, -2 }, 0.5377239,  0.05, -0.03,  -0.3,  -0.6 },
        { { 2,  0,  0, -2, -2, -2 }, 0.5363232,  0.06, -0.03,  -0.4,  -0.7 },
        { { 2, -1,  0, -2,  0, -2 }, 0.5274312,  0.35, -0.20,  -2.4,  -4.1 },
        { { 2,  1,  0, -2, -2, -2 }, 0.5260835,  0.07, -0.04,  -0.5,  -0.8 },
        { { 2,  0,  0, -2,  0, -1 }, 0.5175645, -0.07,  0.04,   0.5,   0.8 },
        { { 2,  0,  0, -2,  0, -2 }, 0.5175251,  1.75, -1.01, -12.2, -21.3 },
        { { 2,  1,  0, -2,  0, -2 }, 0.5079842, -0.05,  0.03,   0.3,   0.6 },
        { { 2,  0, -1, -2,  2, -2 }, 0.5006854,  0.04, -0.03,  -0.3,  -0.6 },
        { { 2,  0,  0, -2,  2, -2 }, 0.5000000,  0.76, -0.44,  -5.5,  -9.6 },
        { { 2,  0,  0,  0,  0,  0 }, 0.4986348,  0.21, -0.12,  -1.5,  -2.6 },
        { { 2,  0,  0,  0,  0, -1 }, 0.4985982,  0.06, -0.04,  -0.4,  -0.8}
    };
    constexpr int M { sizeof(x) / sizeof(x[0]) };
    static_assert( M == 11, "Invalid size for quasi semidiurnal terms in utlibr." );

    //  Compute the harmonic model of dUT1 and dLOD 
    //+ dUT1 and dLOD are set to zero first
    dut1 = dlod = .0e0;

    //  Evaluate the vector of the fundamental arguments
    //+ arg(1:6) = [ GMST+pi, el, elp, f, d, om ] at t = rmjd

    // Convert the input epoch to Julian centuries of TDB since J2000
    double t { (rmjd-RMJD0) / 36525e0 };

    // Compute GMST + pi
    double gmst { std::fmod(67310.54841e0 +
                    t*( (8640184.812866e0 + 3155760000e0) +
                    t*( 0.093104e0 +
                    t*( -0.0000062 ))), 86400e0 ) };

    // Fundamental arguments
    double fargs[6];
    iers2010::fundarg(t, fargs+1);
    fargs[0] = gmst / RAD2SEC + PI;
    fargs[0] = std::fmod(fargs[0], TWOPI);

    double angle, sina, cosa;
    for (int j=0; j<M; j++) {
        //  For the j-th term of the trigonometric expansion, compute the angular
        //+ argument angle of sine and cosine functions as a linear integer
        //+ combination of the 6 fundamental arguments
        angle = .0e0;
        for (int i=0; i<6; i++) {
            angle += ( double(x[j].iarg[i]) * fargs[i] );
        }
        angle = std::fmod(angle, TWOPI);
        // Compute contribution from the j-th term of expansion to dUT1 and dLOD
        sina = std::sin(angle);
        cosa = std::cos(angle);
        dut1 += x[j].dut1s * sina + x[j].dut1c * cosa;
        dlod += x[j].dlods * sina + x[j].dlodc * cosa;
    }

    // Finished.
    return 0;
}
