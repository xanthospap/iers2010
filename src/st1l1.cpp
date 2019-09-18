#include "iers2010.hpp"

inline double
inner_product3(const double* x1, const double* x2) noexcept
{
  return x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
}

/**
 * @details This function gives the corrections induced by the latitude 
 *          dependence given by L^1 in Mathews et al. 1991 (See References).
 *          This function is a translation/wrapper for the fortran ST1L1
 *          subroutine, found here : 
 *          http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  xsta    Geocentric position of the IGS station (Note 1)
 * @param[in]  xsun    Geocentric position of the Sun (Note 2)
 * @param[in]  xmon    Geocentric position of the Moon (Note 2)
 * @param[in]  fac2sun Degree 2 TGP factor for the Sun (Note 3)      
 * @param[in]  fac2mon Degree 2 TGP factor for the Moon (Note 3) 
 * @param[out] xcorsta Out of phase station corrections for semi-diurnal band
 * 
 * @note
 *    -# The IGS station is in ITRF co-rotating frame. All coordinates are
 *       expressed in meters, as arrays, i.e. [x,y,z].
 *    -# The position is in Earth Centered Earth Fixed (ECEF) frame.  All
 *       coordinates are expressed in meters, as arrays, i.e. [x,y,z].
 *    -# The expressions are computed in the main program. TGP is the tide
 *       generated potential. The units are inverse meters.
 *    -# Status: Class 1
 *    -# This fucnction is part of the package dehanttideinel, see
 *       ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/
 *
 * @version 31.07.2009
 * 
 * @cite iers2010,
 *       Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, "Tidal station
 *       displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477,
 *       Mathews, P. M., Buffett, B. A., Herring, T. A., Shapiro, I. I.,
 *       1991b, Forced nutations of the Earth: Influence of inner core
 *       Dynamics 2. Numerical results and comparisons, J. Geophys. Res.,
 *       96, 8243-8257
 * 
 */
void
iers2010::dhtide::st1l1(const double* xsta,const double* xsun,
    const double* xmon, double fac2sun, double fac2mon, double* xcorsta)
{
  const double l1d  { 0.0012e0 },
        l1sd { 0.0024e0 };

  // Compute the normalized position vector of the IGS station.
  const double rsta     { sqrt(inner_product3(xsta,xsta)) };

  const double sinphi   { xsta[2]/rsta };
  const double sinphi2  { sinphi*sinphi };
  const double cosphi   { sqrt(xsta[0]*xsta[0]+xsta[1]*xsta[1])/rsta };
  const double cosphi2  { cosphi*cosphi };
  const double sinla    { xsta[1]/cosphi/rsta };
  const double cosla    { xsta[0]/cosphi/rsta };

  // Compute the normalized position vector of the Moon.
  const double rmon2   { inner_product3(xmon, xmon) };

  // Compute the normalized position vector of the Sun.
  const double rsun2   { inner_product3(xsun, xsun) };

  // Compute the station corrections for the diurnal band.
  double l1    { l1d };
  double dnsun { -l1*sinphi2*fac2sun*xsun[2]*(xsun[0]*cosla+xsun[1]*sinla)
    /rsun2 };
  double dnmon { -l1*sinphi2*fac2mon*xmon[2]*(xmon[0]*cosla+xmon[1]*sinla)
    /rmon2 };
  double desun {  l1*sinphi*(cosphi2-sinphi2)*fac2sun*xsun[2]*
    (xsun[0]*sinla-xsun[1]*cosla)/rsun2 };
  double demon {  l1*sinphi*(cosphi2-sinphi2)*fac2mon*xmon[2]*
    (xmon[0]*sinla-xmon[1]*cosla)/rmon2 };

  double de { 3e0*(desun+demon) };
  double dn { 3e0*(dnsun+dnmon) };

  xcorsta[0] = -de*sinla-dn*sinphi*cosla;
  xcorsta[1] =  de*cosla-dn*sinphi*sinla;
  xcorsta[2] =  dn*cosphi;

  // Compute the station corrections for the semi-diurnal band.
  l1 = l1sd;  
  const double costwola { cosla*cosla-sinla*sinla };
  const double sintwola { 2.e0*cosla*sinla };

  dnsun = -l1/2e0*sinphi*cosphi*fac2sun*((pow(xsun[0],2)-pow(xsun[1],2))*
      costwola+2e0*xsun[0]*xsun[1]*sintwola)/rsun2;

  dnmon = -l1/2e0*sinphi*cosphi*fac2mon*((pow(xmon[0],2)-pow(xmon[1],2))*
      costwola+2e0*xmon[0]*xmon[1]*sintwola)/rmon2;

  desun = -l1/2e0*sinphi2*cosphi*fac2sun*((pow(xsun[0],2)-pow(xsun[1],2))*
      sintwola-2e0*xsun[0]*xsun[1]*costwola)/rsun2;

  demon = -l1/2e0*sinphi2*cosphi*fac2mon*((pow(xmon[0],2)-pow(xmon[1],2))*
      sintwola-2e0*xmon[0]*xmon[1]*costwola)/rmon2;

  de = 3e0*(desun+demon);
  dn = 3e0*(dnsun+dnmon);

  xcorsta[0] += (-de*sinla-dn*sinphi*cosla);
  xcorsta[1] += (de*cosla-dn*sinphi*sinla);
  xcorsta[2] += (dn*cosphi);

  // Finished
  return;
}

/*
#include <stdio.h>
int main ()
{
  double xsta[] = {4075578.385e0,931852.890e0,4801570.154e0};
  double xsun[] = {137859926952.015e0,54228127881.4350e0,23509422341.6960e0};
  double xmon[] = {-179996231.920342e0,-312468450.131567e0,-169288918.592160e0};
  double fac2sun= 0.163271964478954e0;
  double fac2mon= 0.321989090026845e0;
  double xcorsta[3] = {.0e0,.0e0,.0e0};

  st1l1 ( xsta,xsun,xmon,fac2sun,fac2mon,xcorsta );

  printf ("\nSubroutine ST1L1");
  printf ("\nDifferences in meters:");
  printf ("\n\t|dx| = %15.12f", 0.2367189532359759044e-03-xcorsta[0]);
  printf ("\n\t|dy| = %15.12f", 0.5181609907284959182e-03-xcorsta[1]);
  printf ("\n\t|dz| = %15.12f",-0.3014881422940427977e-03-xcorsta[2]);
  printf ("\n");
  
  return 0;
}
*/
