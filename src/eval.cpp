#include "hardisp.hpp"

/**
 * @details This function performs cubic spline interpolation of a given function
 *          sampled at unequally spaced intervals.  The subroutine SPLINE needs
 *          to be called beforehand to set up the array s.
 * 
 * @param[in]  y    The coordinate at which a function value is desired (Note 1)
 * @param[in]  nn   Number of samples of the original function
 * @param[in]  x    Array containing sample coordinates x(1),x(2),...x(nn) (Note 2)
 * @param[in]  s    Array containing the 2nd derivatives at the sample points (Note 3)
 * @param[out] u    Array containing samples of a function at the coordinates 
 *                  x(1),x(2),...x(nn)
 * @param[out] val  The interpolated value of the function at y
 * @return          Zero on success; else 1.
 *
 * @note
 *     -# If y falls outside the range (x(1),x(nn)), the value at the nearest
 *        endpoint of the series is used.
 *     -# The sequence x(1),x(2),...x(nn) must be strictly increasing.
 *     -# This array is found by the function SPLINE, which must be called
 *        once before beginning this interpolation.
 * 
 * @version 2009 August 19
 * 
 */
int iers2010::hisp::eval (const double& y,const int& nn,const double* x,
        const double* u,const double* s,double& val)
{
  if (nn < 0)
      return 1;
  
  // find index k, such that x[k-1] < y <= x[k]
  int k = std::distance (x,std::lower_bound (x,x+nn,y));
  
  int k2 = k;
  int k1 = k-1;
  
  // handle out of bounds; substitute endpoint values
  if (k >= nn) {
    //double dk = x[k2] - x[k1];
    val = u[nn-1];
    return 0;
  } 
  else if (k == 0) {
    val = u[0];
    return 0;
  }
  
  //  Evaluate and then interpolate.
  //+ Note that this can fail if dk is ~0
  double dy = x[k2] - y;
  double dy1 = y - x[k1];
  double dk = x[k2] - x[k1];
  if (dk < 1e-30)
      return 1;
  double deli = 1.0e0 / (6.0e0 * dk);
  double ff1 = s[k1]*dy*dy*dy;
  double ff2 = s[k2]*dy1*dy1*dy1;
  double f1  = (ff1+ff2) * deli;
  double f2 = dy1*((u[k2]/dk)-(s[k2]*dk)/6.0e0);
  double f3 = dy*((u[k1]/dk)-(s[k1]*dk)/6.0e0);
  val = f1 + f2 + f3;
  //printf ("\nEVAL: F1=%14.6f F2=%14.6f F3=%14.6f",f1,f2,f3);
  
  // Finished
  return 0;
}