// https://vmf.geo.tuwien.ac.at/codes/saasthyd.m
#include "tropo.hpp"
#include <cmath>

// This subroutine determines the zenith hydrostatic delay based on the
// equation by Saastamoinen (1972) as refined by Davis et al. (1985)
//
// c Reference:
// Saastamoinen, J., Atmospheric correction for the troposphere and
// stratosphere in radio ranging of satellites. The use of artificial
// satellites for geodesy, Geophys. Monogr. Ser. 15, Amer. Geophys. Union,
// pp. 274-251, 1972.
// Davis, J.L, T.A. Herring, I.I. Shapiro, A.E.E. Rogers, and G. Elgered,
// Geodesy by Radio Interferometry: Effects of Atmospheric Modeling Errors
// on Estimates of Baseline Length, Radio Science, Vol. 20, No. 6,
// pp. 1593-1607, 1985.
//
// input parameters:
// p:     pressure in hPa
// dlat:  ellipsoidal latitude in radians
// dlon:  longitude in radians
// hell:  ellipsoidal height in m
//
// output parameters:
// zhd:  zenith hydrostatic delay in meter
//
// Example 1 :
// p = 1000;
// dlat = 48d0*pi/180.d0
// hell = 200.d0
//
// output:
// zhd = 2.2695 m
double dso::saasthyd(double p, double dlat, double hell) noexcept {
  // calculate denominator f
  const double f = 1e0 - 0.00266e0 * std::cos(2e0 * dlat) - 0.00000028e0 * hell;

  // calculate the zenith hydrostatic delay zhd
  return 0.0022768 * p / f;
}
