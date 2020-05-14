#ifndef __IERS_1010_DEHANTTIDEINEL__
#define __IERS_1010_DEHANTTIDEINEL__

#include <numeric>
#include <cmath>
#ifdef DEBUG
#include <iostream>
#include <cfenv>
#endif

namespace iers2010
{
  
namespace dhtide
{
  /// Out-of-phase corrections induced by mantle anelasticity in the diurnal
  /// band.
  void
  st1idiu(const double*, const double*, const double*, double, double, double*)
  noexcept;

  /// Out-of-phase corrections induced by mantle anelasticity in the
  /// semi-diurnal band. 
  void
  st1isem(const double*, const double*, const double*, double, double, double*)
  noexcept;

  /// Corrections induced by the latitude dependence given by L^1 in 
  /// Mathews et al. 1991. 
  void
  st1l1(const double*, const double*, const double*, double, double, double*)
  noexcept;

  /// In-phase and out-of-phase corrections induced by mantle anelasticity
  /// in the diurnal band.
  void
  step2diu(const double*, double, double, double*) noexcept;

  /// In-phase and out-of-phase corrections induced by mantle anelasticity
  /// in the long period band.
  void
  step2lon(const double*, double, double*) noexcept;

  /// @brief Function to compute the scalar product of two vectors and 
  ///        their norms.
  /// @note  If the vectors contain more than 3 elements, only he first 3 are
  ///        used for the computation.
  /// @param[in]  x  Vector of dimension (at least) 3
  /// @param[in]  y  Vector of dimension (at least) 3
  /// @param[out] r1 (Euclidean) norm of vector x
  /// @param[out] r2 (Euclidean) norm of vector y
  /// @return        scalar product of vectors x and y
  inline double
  sprod(const double* x, const double* y, double& r1, double& r2)
  noexcept
  {
    r1 = std::sqrt(std::inner_product(x,x+3,x,.0e0));
    r2 = std::sqrt(std::inner_product(y,y+3,y,.0e0));
    return std::inner_product(x, x+3, y, .0e0);
  }

}// dhtide
  
}// iers2010
#endif
