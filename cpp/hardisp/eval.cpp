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
 * @param[in]  u    Array containing samples of a function at the coordinates 
 *                  x(1),x(2),...x(nn)
 * @param[out] eval The interpolated value of the function at y
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
int iers2010::hisp::eval (const double& y, const int& nn, const double* x, const double* s, const double* u, double& eval)
{
  if (nn < 0) return 1;
    
  // If y is out of range, substitute endpoint values
  if (y <= x[0]) {
    eval = x[0];
    return 0;
  }
  if (y >= x[nn-1]) {
    eval = u[nn-1];
    return 0;
  }
  
  //  Locate interval (x(k1),x(k2)) which contains y
  //+ Keep it simple; use the std algorithm (binary search)
  //+ since the array is already sorted in ascending order.
  DO 100 K=2,NN
  IF(X(K-1).LT.Y.AND.X(K).GE.Y) THEN
  K1=K-1
  K2=K
  ENDIF
  100   CONTINUE
}