#include <numeric>
#include "iers2010.hpp"

/**
 * @details This function gives the in-phase and out-of-phase corrections
 *          induced by mantle anelasticity in the diurnal band.
 *          This function is a translation/wrapper for the fortran STEP2DIU
 *          subroutine, found here at 
 *          http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  xsta    Geocentric position of the IGS station (Note 1)
 * @param[in]  fhr     Fractional hours in the day (Note 2)
 * @param[in]  t       Centuries since J2000
 * @param[out] xcorsta In phase and out of phase station corrections
 *                     for diurnal band (Note 4) 
 * 
 * @note
 *    -# The IGS station is in ITRF co-rotating frame. All coordinates are
 *       expressed in meters, as arrays, i.e. [x,y,z].
 *    -# The fractional hours in the day is computed as the hour + minutes/60.0
 *       + sec/3600.0.  The unit is expressed in Universal Time (UT).
 *    -# ----
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
iers2010::dhtide::step2diu(const double* xsta, double fhr, double t,
    double* xcorsta)
{

  // Set constants
#ifdef USE_EXTERNAL_CONSTS
#else
  constexpr double D2PI   ( 6.283185307179586476925287e0 );  // 2*pi
  constexpr double DEG2RAD( D2PI / 360e0 );   // degrees to radians
#endif

  constexpr double datdi[][9] = {
    {-3e0,  0e0,  2e0,  0e0,  0e0, -0.01e0,     0e0,    0e0,    0e0},   
    {-3e0,  2e0,  0e0,  0e0,  0e0, -0.01e0,     0e0,    0e0,    0e0},   
    {-2e0,  0e0,  1e0, -1e0,  0e0, -0.02e0,     0e0,    0e0,    0e0},   
    {-2e0,  0e0,  1e0,  0e0,  0e0, -0.08e0,     0e0,-0.01e0, 0.01e0},
    {-2e0,  2e0, -1e0,  0e0,  0e0, -0.02e0,     0e0,    0e0,    0e0},
    {-1e0,  0e0,  0e0, -1e0,  0e0, -0.10e0,     0e0,    0e0,    0e0},
    {-1e0,  0e0,  0e0,  0e0,  0e0, -0.51e0,     0e0,-0.02e0, 0.03e0},
    {-1e0,  2e0,  0e0,  0e0,  0e0,  0.01e0,     0e0,    0e0,    0e0},
    { 0e0, -2e0,  1e0,  0e0,  0e0,  0.01e0,     0e0,    0e0,    0e0},
    { 0e0,  0e0, -1e0,  0e0,  0e0,  0.02e0,     0e0,    0e0,    0e0},
    { 0e0,  0e0,  1e0,  0e0,  0e0,  0.06e0,     0e0,    0e0,    0e0},
    { 0e0,  0e0,  1e0,  1e0,  0e0,  0.01e0,     0e0,    0e0,    0e0},
    { 0e0,  2e0, -1e0,  0e0,  0e0,  0.01e0,     0e0,    0e0,    0e0},
    { 1e0, -3e0,  0e0,  0e0,  1e0, -0.06e0,     0e0,    0e0,    0e0},
    { 1e0, -2e0,  0e0, -1e0,  0e0,  0.01e0,     0e0,    0e0,    0e0},
    { 1e0, -2e0,  0e0,  0e0,  0e0, -1.23e0, -0.07e0, 0.06e0, 0.01e0},
    { 1e0, -1e0,  0e0,  0e0, -1e0,  0.02e0,     0e0,    0e0,    0e0},
    { 1e0, -1e0,  0e0,  0e0,  1e0,  0.04e0,     0e0,    0e0,    0e0},
    { 1e0,  0e0,  0e0, -1e0,  0e0, -0.22e0,  0.01e0, 0.01e0,    0e0},
    { 1e0,  0e0,  0e0,  0e0,  0e0, 12.00e0, -0.80e0,-0.67e0,-0.03e0},
    { 1e0,  0e0,  0e0,  1e0,  0e0,  1.73e0, -0.12e0,-0.10e0,    0e0},
    { 1e0,  0e0,  0e0,  2e0,  0e0, -0.04e0,     0e0,    0e0,    0e0}, 
    { 1e0,  1e0,  0e0,  0e0, -1e0, -0.50e0, -0.01e0, 0.03e0,    0e0},
    { 1e0,  1e0,  0e0,  0e0,  1e0,  0.01e0,     0e0,    0e0,    0e0},
    { 0e0,  1e0,  0e0,  1e0, -1e0, -0.01e0,     0e0,    0e0,    0e0},
    { 1e0,  2e0, -2e0,  0e0,  0e0, -0.01e0,     0e0,    0e0,    0e0},
    { 1e0,  2e0,  0e0,  0e0,  0e0, -0.11e0,  0.01e0, 0.01e0,    0e0},
    { 2e0, -2e0,  1e0,  0e0,  0e0, -0.01e0,     0e0,    0e0,    0e0},
    { 2e0,  0e0, -1e0,  0e0,  0e0, -0.02e0,     0e0,    0e0,    0e0},
    { 3e0,  0e0,  0e0,  0e0,  0e0,     0e0,     0e0,    0e0,    0e0},
    { 3e0,  0e0,  0e0,  1e0,  0e0,     0e0,     0e0,    0e0,    0e0}
  };

  // Compute the phase angles in degrees.
  double s  {   218.31664563e0
          + (481267.88194e0
          + (    -0.0014663889e0 
          + (     0.00000185139e0)*t)*t)*t 
          };

  double tau { fhr*15e0
          +    280.4606184e0
          + (36000.7700536e0
          + (    0.00038793e0
          + (   -0.0000000258e0)*t)*t)*t
          + (-s)
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
          + (    0.00000001650e0)*t)*t)*t)*t
          };

  double ps { 282.93734098e0
          + (   1.71945766667e0
          + (   0.00045688889e0
          + (  -0.00000001778e0
          + (  -0.00000000334e0)*t)*t)*t)*t
          };

  // Reduce angles to between the range 0 and 360.
  s   = std::fmod(s,   360e0);
  tau = std::fmod(tau, 360e0);
  h   = std::fmod(h,   360e0);
  p   = std::fmod(p,   360e0);
  zns = std::fmod(zns, 360e0);
  ps  = std::fmod(ps,  360e0);

  const double rsta   { std::sqrt(std::inner_product(
        xsta, xsta+3, xsta, 0e0)) };
  const double sinphi { xsta[2]/rsta };
  const double cosphi { std::sqrt(xsta[0]*xsta[0]+xsta[1]*xsta[1])/rsta };

  const double cosla  { xsta[0]/cosphi/rsta };
  const double sinla  { xsta[1]/cosphi/rsta };
  const double zla    { std::atan2(xsta[1], xsta[0]) };

  xcorsta[0] = xcorsta[1] = xcorsta[2] = 0e0;

  double thetaf, dr, dn, de;
  for (int j=0; j<31; j++) {
    // Convert from degrees to radians.
    thetaf = (tau+datdi[j][0]*s+datdi[j][1]*h+datdi[j][2]*p+
        datdi[j][3]*zns+datdi[j][4]*ps)*DEG2RAD;

    dr     = datdi[j][5]*2e0*sinphi*cosphi*sin(thetaf+zla)+
      datdi[j][6]*2e0*sinphi*cosphi*cos(thetaf+zla);

    dn     = datdi[j][7]*(cosphi*cosphi-sinphi*sinphi)*sin(thetaf+zla)+
      datdi[j][8]*(cosphi*cosphi-sinphi*sinphi)*cos(thetaf+zla);

    // DE=DATDI(8,J)*SINPHI*COS(THETAF+ZLA)+
    // Modified 20 June 2007
    de     =  datdi[j][7]*sinphi*cos(thetaf+zla)-
      datdi[j][8]*sinphi*sin(thetaf+zla);

    xcorsta[0] += dr*cosla*cosphi-de*sinla-dn*sinphi*cosla;  
    xcorsta[1] += dr*sinla*cosphi+de*cosla-dn*sinphi*sinla;
    xcorsta[2] += +dr*sinphi+dn*cosphi;
  }   

  xcorsta[0] /= 1000e0;
  xcorsta[1] /= 1000e0;
  xcorsta[2] /= 1000e0;

  // Finished
  return;
}

/*
#include <stdio.h>
int main ()
{
  double xsta[] = {4075578.385e0,931852.890e0,4801570.154e0};
  double fhr    = .0e0;
  double t      = 0.1059411362080767e0;
  double xcorsta[3] = {.0e0,.0e0,.0e0};

  step2diu ( xsta,fhr,t,xcorsta );

  printf ("\nSubroutine STEP2DIU");
  printf ("\nDifferences in meters:");
  printf ("\n\t|dx| = %15.12f",0.4193085327321284701e-02 -xcorsta[0]);
  printf ("\n\t|dy| = %15.12f",0.1456681241014607395e-02 -xcorsta[1]);
  printf ("\n\t|dz| = %15.12f",0.5123366597450316508e-02 -xcorsta[2]);
  printf ("\n");
  
  return 0;
}
*/
