#include <numeric>
#include "iers2010.hpp"

/**
 * @details This function gives the in-phase and out-of-phase corrections
 *          induced by mantle anelasticity in the long period band.
 *          This function is a translation/wrapper for the fortran STEP2LON
 *          subroutine, found here : 
 *          http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  xsta    Geocentric position of the IGS station (Note 1)
 * @param[in]  t       Centuries since J2000
 * @param[out] xcorsta In phase and out of phase station corrections
 *                     for diurnal band (Note 2)
 * 
 * @note
 *    -# The IGS station is in ITRF co-rotating frame. All coordinates are
 *       expressed in meters, as arrays, i.e. [x,y,z].
 *    -# All coordinates are expressed in meters.
 *    -# Status: Class 1
 *    -# This fucnction is part of the package dehanttideinel, see
 *       ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/
 *
 * @version 20.10.2010
 * 
 * @cite iers2010,
 *       Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, "Tidal station
 *       displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477,
 * 
 */
void
iers2010::dhtide::step2lon(const double* xsta, double t, double* xcorsta)
{
    // Set constants
    #ifdef USE_EXTERNAL_CONSTS
        /*constexpr double TWOPI   (D2PI);*/
    #else
        constexpr double D2PI   ( 6.283185307179586476925287e0 );// 2*pi
        constexpr double DEG2RAD( D2PI / 360e0 );                // degrees to radians
    #endif

    /*static*/ const double datdi[][9] = {
        {0e0, 0e0, 0e0, 1e0, 0e0,   0.47e0,  0.23e0,  0.16e0,  0.07e0},
        {0e0, 2e0, 0e0, 0e0, 0e0,  -0.20e0, -0.12e0, -0.11e0, -0.05e0},
        {1e0, 0e0,-1e0, 0e0, 0e0,  -0.11e0, -0.08e0, -0.09e0, -0.04e0},
        {2e0, 0e0, 0e0, 0e0, 0e0,  -0.13e0, -0.11e0, -0.15e0, -0.07e0},
        {2e0, 0e0, 0e0, 1e0, 0e0,  -0.05e0, -0.05e0, -0.06e0, -0.03e0}
    };

    // Compute the phase angles in degrees.
    double s {    218.31664563e0
            + (481267.88194e0
            + (    -0.0014663889e0 
            + (     0.00000185139e0)*t)*t)*t
            }; 

    double pr { (1.396971278e0
            + (  0.000308889e0
            + (  0.000000021e0
            + (  0.000000007e0)*t)*t)*t)*t
            };

    s += pr;

    double h {   280.46645e0
            + (36000.7697489e0
            + (    0.00030322222e0 
            + (    0.000000020e0
            + (   -0.00000000654e0)*t)*t)*t)*t
            };

    double p {   83.35324312e0
            + (4069.01363525e0
            + (  -0.01032172222e0
            + (  -0.0000124991e0
            + (   0.00000005263e0)*t)*t)*t)*t
            };

    double zns { 234.95544499e0
            + ( 1934.13626197e0
            + (   -0.00207561111e0
            + (   -0.00000213944e0
            + (   0.00000001650e0)*t)*t)*t)*t
            };

    double ps { 282.93734098e0
            + (   1.71945766667e0
            + (   0.00045688889e0
            + (  -0.00000001778e0
            + (  -0.00000000334e0)*t)*t)*t)*t
            };

    // Reduce angles to between the range 0 and 360.
    s   = std::fmod(s,   360e0);
    h   = std::fmod(h,   360e0);
    p   = std::fmod(p,   360e0);
    zns = std::fmod(zns, 360e0);
    ps  = std::fmod(ps,  360e0);

    const double rsta   { std::sqrt(std::inner_product(
        xsta, xsta+3, xsta, .0e0)) };
    const double sinphi { xsta[2]/rsta };
    const double cosphi { std::sqrt(xsta[0]*xsta[0]+xsta[1]*xsta[1])/rsta };

    const double cosla  { xsta[0]/cosphi/rsta };
    const double sinla  { xsta[1]/cosphi/rsta };

    double dr_tot { .0e0 },
           dn_tot { .0e0 };

    xcorsta[0] = xcorsta[1] = xcorsta[2] = .0e0;

    double thetaf, dr, dn, de;
    for (int j = 0; j < 5; j++) {

        thetaf = (datdi[j][0]*s+datdi[j][1]*h+datdi[j][2]*p+
                datdi[j][3]*zns+datdi[j][4]*ps)*DEG2RAD;

        dr     = datdi[j][5]*(3e0*sinphi*sinphi-1e0)/2e0*cos(thetaf)+
                datdi[j][7]*(3e0*sinphi*sinphi-1e0)/2e0*sin(thetaf);

        dn     = datdi[j][6]*(cosphi*sinphi*2e0)*cos(thetaf)+
                datdi[j][8]*(cosphi*sinphi*2e0)*sin(thetaf);

        de     = 0e0;
 
        dr_tot += dr;
        dn_tot += dn;

        xcorsta[0] += dr*cosla*cosphi-de*sinla-dn*sinphi*cosla;
        xcorsta[1] += dr*sinla*cosphi+de*cosla-dn*sinphi*sinla;  
        xcorsta[2] += dr*sinphi+dn*cosphi;
    }

    xcorsta[0] /= 1000e0;
    xcorsta[1] /= 1000e0;
    xcorsta[2] /= 1000e0;

    // Finished.
    return;
}

/*
#include <stdio.h>
int main ()
{
  double xsta[] = {4075578.385e0,931852.890e0,4801570.154e0};
  double t      = 0.1059411362080767e0;
  double xcorsta[3] = {.0e0,.0e0,.0e0};

  step2lon ( xsta,t,xcorsta );

  printf ("\nSubroutine STEP2LON");
  printf ("\nDifferences in meters:");
  printf ("\n\t|dx| = %15.12f",-0.9780962849562107762e-04 -xcorsta[0]);
  printf ("\n\t|dy| = %15.12f",-0.2236349699932734273e-04 -xcorsta[1]);
  printf ("\n\t|dz| = %15.12f",0.3561945821351565926e-03 -xcorsta[2]);
  printf ("\n");
  
  return 0;
}
*/
