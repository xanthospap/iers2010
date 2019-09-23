#ifndef __NGPT_IERS_1010__
#define __NGPT_IERS_1010__
#include <cmath>

namespace iers2010
{
  /// Compute the lunisolar fundamental arguments.
  int
  fundarg(double, double*);

  /// Compute the diurnal lunisolar effect on polar motion.
  int
  pmsdnut2(double, double&, double&);

  /// Compute the subdiurnal librations in UT1.
  int
  utlibr(double, double&, double&);

  /// Compute corrections to the coordinates of the CIP to account for 
  /// Free Core Nutation.
  int
  fcnnut(double, double&, double&, double&, double&);

  /// Compute the angular argument which depends on time for 11 tidal
  /// argument calculations.
  int
  arg2(int, double, double*);

  /// Compute tidal corrections of station displacements caused by lunar and
  /// solar gravitational attraction.
  int
  dehanttideinel(const double* xsta,const double* xsun, const double* xmon,
      int yr, int month, int day, double fhr, double* dxtide);

  /// Compute the diurnal and semi-diurnal variations in Earth Orientation
  /// Parameters from ocean tides.
  int
  ortho_eop(double, double&, double&, double&);

  /// Evaluate the effects of zonal Earth tides on the rotation of the Earth.
  int
  rg_zont2(double, double&, double&, double&);

  /// Compute the global total FCULa mapping function.
  double
  fcul_a(double, double, double, double);

  /// Computes the global total FCULb mapping function.
  double
  fcul_b(double, double, double, double);

  // Determine the total zenith delay following Mendes and Pavlis, 2004. 
  int
  fcul_zd_hpa(double, double, double, double, double, double&, double&,
    double&);

  /// Determines the total zenith delay following (Mendes and Pavlis, 2004).
  int
  fculzd_hpa(double, double, double, double, double, double&, double&,
    double&);

  /// Compute the Global Mapping Functions (GMF).
  int
  gmf(double, double, double, double, double, double&, double&);

  /// Compute the Vienna Mapping Functions 1 (VMF1), to be used with "a" 
  /// coefficients computed for a given site.
  int
  vmf1(double, double, double, double, double, double&, double&);

  /// Compute the Vienna Mapping Functions 1 (VMF1), with height corrections,
  /// to be used with "a" coefficients computed for a grid.
  int
  vmf1_ht(double, double, double, double, double, double, double&, double&);

  /// Compute the Global Pressure and Temperature (GPT), based on spherical 
  /// harmonics up to degree and order 9. 
  int
  gpt(double, double, double, double, double&, double&, double&);

  /// Compute the Global Pressure and Temperature 2 model (GPT2), 
  /// and the "a" coefficients for vmf1_ht.
  int
  gpt2(double, double*, double*, double*, int, int, double*, double*, 
    double*, double*, double*, double*, double*, const char* ifile=nullptr);
  
  namespace oeop
  {
      /// Compute the time dependent part of second degree diurnal and 
      /// semidiurnal tidal potential.
      int
      cnmtx(double dmjd, double* h);

  } // oeop

} // iers2010

#endif
