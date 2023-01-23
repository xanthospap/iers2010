#ifndef __CUBIC_SPLINE_INTERPOLATION_DSO_HPP__
#define __CUBIC_SPLINE_INTERPOLATION_DSO_HPP__

#include <algorithm>
#include <cassert>

namespace dso {

/// @brief Cubic Spline interpolation. This function assumes the second
///        derivative of the interpolating function has already been computed
///        for the data points (x,y) and stored in the y2 array (see function
///        cspline_deriv).
///
///        We want to compute the interpolated value at atx, given the
///        data points stored in the x and y arrays (of size n) and the
///        seconds derivatives (alreadt computed) available at the size-n
///        array y2. The interpolated y-value is stored at yintrp, if and only
///        if the return value is 0.
///        If (at input) we do not know the index atx has in the input x-data
///        array, set index to < 0. Else, the function will find the index by
///        searching through the x array
///
///        Note that if n < 4, then we will perform linear interpolation.
///        Also, if the interpoated x is out of range, we will return the
///        boundary y values (either y[0] or y[n-1])
///
/// Reference: Numerical Recipes, 3rd edition
///
/// @param[in] atx The x-value we want to interpolate at
/// @param[in] index The index atx value corresponds to, aka:
///            x[index] <= atx < x[index+1]
///            If not known, set to any value less than 0 and the function
///            will find it (at set it)
/// @param[in] x Array of size n, containing x-values
/// @param[in] y Array of size n, containing y-values (corresponding to the
///              array above)
/// @param[in] n Size of x and y input arrays, and also of y2 array
/// @param[in] y2 Array of size n, containing derivatives of interpolating
///              function at given nodes
/// @param[out] yintrp The interpolated y value at atx
/// @return Anything greater than 0 denotes an error. A value of -1 denotes
///         that no interpolation is performed, but boundary values are 
///         extrapolated. If -2 is returned, the function has performed linear
///         interpolation. 0 denotes cubic spline interpolation
int cspline_interp(double atx, int &index, const double *const x,
                   const double *const y, int n, const double *const y2,
                   double &yintrp) noexcept;

/// @brief Compute y2[0..n-1] array containing second derivatives of the
/// interpolating function for cuqbic spline interpolation.
/// If yp1 and/or ypn are equal to 1x10^99 or larger, the routine is signaled
/// to set the corresponding boundary condition for a natural spline, with
/// zero second derivative on that boundary; otherwise, they are the values
/// of the ﬁrst derivatives at the endpoints.
///
/// Reference: Numerical Recipes, 3rd edition
///
/// @param[in] x Array of size n, containing x-values
/// @param[in] y Array of size n, containing y-values (corresponding to the
///              a array above)
/// @param[in] n Size of x and y input arrays, and also of y2 output array and
///              u (workspace array). User is responsible for allocating space
///              for u, which should be at least n
/// @param[out] y2 Array of size n, on exit contains derivatives of
///              interpoalting function at given nodes
/// @param[in] yp1 Signal value for value of first derivative at starting
///              point (see description)
/// @param[in] ypn Signal value for value of first derivative at ending
///              point (see description)
void cspline_deriv(const double *const x, const double *const y, int n,
                   double *__restrict__ y2, double yp1, double ypn,
                   double *__restrict__ u) noexcept;
} // namespace dso

#endif
