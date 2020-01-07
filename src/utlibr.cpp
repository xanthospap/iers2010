#include "iers2010.hpp"
#ifdef USE_EXTERNAL_CONSTS
#include "gencon.hpp"
#endif

/// @details This function evaluates the model of subdiurnal libration
///          in the axial component of rotation, expressed by UT1 and LOD.
///          This effect is due to the influence of tidal gravitation on the
///          departures of the Earth's mass distribution from the rotational
///          symmetry, expressed by the non-zonal components of geopotential.
///          The amplitudes have been computed for an elastic Earth with liquid
///          core. The adopted truncation level is 0.033 microseconds in UT1
///          corresponding to the angular displacement of 0.5 microarcseconds
///          or to 0.015 mm at the planet surface. With this truncation level
///          the model contains 11 semidiurnal terms. The coefficients of
///          the model are given in Table 5.1b of the IERS Conventions (2010).
///          This function is a translation/wrapper for the fortran UTLIBR
///          subroutine, found here : 
///          http://maia.usno.navy.mil/conv2010/software.html
/// 
/// @param[in]  rmjd Time expressed as modified Julian date
/// @param[out] dut1 Incremental UT1 in microseconds
/// @param[out] dlod Incremental LOD in microseconds per day
/// @return          An integer value, always 0. 
///
/// @note
///   - The procedure fundarg is the same as used by the program pmsdnut2
///   which implements the corresponding model of the lunisolar libration in
///   polar motion.
///   - Calls function iers2010::fundarg
/// 
/// @version 23.06.2010
/// 
/// @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS 
///       Technical Note No. 36, BKG (2010); Chapter 5.5.3.3
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
    {{ 2,-2, 0,-2, 0,-2 }, 0.5377239e0, 0.05e0,-0.03e0, -0.3e0, -0.6e0},
    {{ 2, 0, 0,-2,-2,-2 }, 0.5363232e0, 0.06e0,-0.03e0, -0.4e0, -0.7e0},
    {{ 2,-1, 0,-2, 0,-2 }, 0.5274312e0, 0.35e0,-0.20e0, -2.4e0, -4.1e0},
    {{ 2, 1, 0,-2,-2,-2 }, 0.5260835e0, 0.07e0,-0.04e0, -0.5e0, -0.8e0},
    {{ 2, 0, 0,-2, 0,-1 }, 0.5175645e0,-0.07e0, 0.04e0,  0.5e0,  0.8e0},
    {{ 2, 0, 0,-2, 0,-2 }, 0.5175251e0, 1.75e0,-1.01e0,-12.2e0,-21.3e0},
    {{ 2, 1, 0,-2, 0,-2 }, 0.5079842e0,-0.05e0, 0.03e0,  0.3e0,  0.6e0},
    {{ 2, 0,-1,-2, 2,-2 }, 0.5006854e0, 0.04e0,-0.03e0, -0.3e0, -0.6e0},
    {{ 2, 0, 0,-2, 2,-2 }, 0.5000000e0, 0.76e0,-0.44e0, -5.5e0, -9.6e0},
    {{ 2, 0, 0, 0, 0, 0 }, 0.4986348e0, 0.21e0,-0.12e0, -1.5e0, -2.6e0},
    {{ 2, 0, 0, 0, 0,-1 }, 0.4985982e0, 0.06e0,-0.04e0, -0.4e0, -0.8e0}
  };
  constexpr int M { sizeof(x) / sizeof(x[0]) };
  static_assert( M == 11, "Invalid size for quasi semidiurnal terms in utlibr." );

  //  Compute the harmonic model of dUT1 and dLOD 
  //+ dUT1 and dLOD are set to zero first
  dut1 = dlod = .0e0;

  //  Evaluate the vector of the fundamental arguments
  //+ arg(1:6) = [ GMST+pi, el, elp, f, d, om ] at t = rmjd

  // Convert the input epoch to Julian centuries of TDB since J2000
  const double t { (rmjd-RMJD0) / 36525e0 };

  // Compute GMST + pi
  const double gmst { std::fmod(67310.54841e0 +
                  t*((8640184.812866e0 + 3155760000e0) +
                  t*(0.093104e0 +
                  t*(-0.0000062e0))), 86400e0) };

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
    angle = 0e0;
    for (int i=0; i<6; i++) angle += (double(x[j].iarg[i]) * fargs[i]);
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
