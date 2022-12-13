#include "datetime/dtfund.hpp"
#include "geodesy/units.hpp"
#include "iers2010.hpp"
#include <datetime/dtcalendar.hpp>

namespace {
// Coefficients of the long and quasi diurnal periodic terms in polar motion
constexpr const struct {
  int iarg[6];
  double per, xs, xc, ys, yc;
} x[] = {
    //  Coefficients of the long periodic terms in polar motion
    //+ Source: IERS Conventions (2010), Table 5.1a
    {{0, 0, 0, 0, 0, -1}, 6798.3837e0, 0.0e0, 0.6e0, -0.1e0, -0.1e0},
    {{0, -1, 0, 1, 0, 2}, 6159.1355e0, 1.5e0, 0.0e0, -0.2e0, 0.1e0},
    {{0, -1, 0, 1, 0, 1}, 3231.4956e0, -28.5e0, -0.2e0, 3.4e0, -3.9e0},
    {{0, -1, 0, 1, 0, 0}, 2190.3501e0, -4.7e0, -0.1e0, 0.6e0, -0.9e0},
    {{0, 1, 1, -1, 0, 0}, 438.35990e0, -0.7e0, 0.2e0, -0.2e0, -0.7e0},
    {{0, 1, 1, -1, 0, -1}, 411.80661e0, 1.0e0, 0.3e0, -0.3e0, 1.0e0},
    {{0, 0, 0, 1, -1, 1}, 365.24219e0, 1.2e0, 0.2e0, -0.2e0, 1.4e0},
    {{0, 1, 0, 1, -2, 1}, 193.55971e0, 1.3e0, 0.4e0, -0.2e0, 2.9e0},
    {{0, 0, 0, 1, 0, 2}, 27.431826e0, -0.1e0, -0.2e0, 0.0e0, -1.7e0},
    {{0, 0, 0, 1, 0, 1}, 27.321582e0, 0.9e0, 4.0e0, -0.1e0, 32.4e0},
    {{0, 0, 0, 1, 0, 0}, 27.212221e0, 0.1e0, 0.6e0, 0.0e0, 5.1e0},
    {{0, -1, 0, 1, 2, 1}, 14.698136e0, 0.0e0, 0.1e0, 0.0e0, 0.6e0},
    {{0, 1, 0, 1, 0, 1}, 13.718786e0, -0.1e0, 0.3e0, 0.0e0, 2.7e0},
    {{0, 0, 0, 3, 0, 3}, 9.1071941e0, -0.1e0, 0.1e0, 0.0e0, 0.9e0},
    {{0, 0, 0, 3, 0, 2}, 9.0950103e0, -0.1e0, 0.1e0, 0.0e0, 0.6e0},
    //  Coefficients of the quasi diurnal terms in polar motion
    //+ Source: IERS Conventions (2010), Table 5.1a
    {{1, -1, 0, -2, 0, -1}, 1.1196992e0, -0.4e0, 0.3e0, -0.3e0, -0.4e0},
    {{1, -1, 0, -2, 0, -2}, 1.1195149e0, -2.3e0, 1.3e0, -1.3e0, -2.3e0},
    {{1, 1, 0, -2, -2, -2}, 1.1134606e0, -0.4e0, 0.3e0, -0.3e0, -0.4e0},
    {{1, 0, 0, -2, 0, -1}, 1.0759762e0, -2.1e0, 1.2e0, -1.2e0, -2.1e0},
    {{1, 0, 0, -2, 0, -2}, 1.0758059e0, -11.4e0, 6.5e0, -6.5e0, -11.4e0},
    {{1, -1, 0, 0, 0, 0}, 1.0347187e0, 0.8e0, -0.5e0, 0.5e0, 0.8e0},
    {{1, 0, 0, -2, 2, -2}, 1.0027454e0, -4.8e0, 2.7e0, -2.7e0, -4.8e0},
    {{1, 0, 0, 0, 0, 0}, 0.9972696e0, 14.3e0, -8.2e0, 8.2e0, 14.3e0},
    {{1, 0, 0, 0, 0, -1}, 0.9971233e0, 1.9e0, -1.1e0, 1.1e0, 1.9e0},
    {{1, 1, 0, 0, 0, 0}, 0.9624365e0, 0.8e0, -0.4e0, 0.4e0, 0.8e0}};

constexpr const int M{sizeof(x) / sizeof(x[0])};
static_assert(M == 25, "Invalid quasi diurnal terms in pmsdnut2.");
}// unnamed namespace

// fargs should have been computed using the compute_fargs function (size=6)
int iers2010::utils::pmsdnut2(const dso::TwoPartDate &mjd,
                              const double *const fargs, double &dx,
                              double &dy) noexcept {
  /*
   *         ----------------------------
   *           D E F I N I T I O N S
   *         ----------------------------
   *
   *  iband  - parameter defining the range of periods for the terms which
   *           are included in computations; if equal to 1 only the quasi
   *           diurnal terms are computed, otherwise the full model
   *  iarg   - array defining for each of the 25 trigonometric terms a set
   *           of 6 integer multipliers of the fundamental angular arguments
   *  arg    - vector of the following 6 fundamental arguments used to
   *           compute the angular argument of the trigonometric functions
   *           arg(1:6) = [ GMST+pi, el, elp, f, d, om ]; this vector is
   *           evaluated by the subroutine FUNDARG which is called as an
   *           external subroutine.  Originally evaluated by the subroutine
   *           PMARGS.
   *  period - array of periods of the trigonometric terms of expansion, in
   *           mean solar days; only for a check - not used in computations
   *  xs, xc - sine and cosine coefficients of the x coordinate of the pole,
   *           in microarcseconds
   *  ys, yc - sine and cosine coefficients of the y coordinate of the pole,
   *           in microarcseconds
   *  angle  - angular argument of the trigonometric functions
   *           angle = Sum(i=1:6) iarg(i,j)*arg(i), for j=1,25
   *
   */
  constexpr const int iband = 1;

  // Rate of secular polar motion, in microarcseconds per year
  // Source: IERS Conventions (2010), Table 5.1a
  constexpr double xrate = -3.8e0, yrate = -4.3e0;

  // Compute the periodical part of the model
  // Coordinates of the pole are set to zero first
  dx = dy = 0e0;

  const int jstart = (iband == 1) ? 15 : 0;
  for (int j = jstart; j < M; j++) {
    // For the j-th term of the trigonometric expansion, compute the angular
    // argument angle of sine and cosine functions as a linear integer
    // combination of the 6 fundamental arguments
    double angle = 0e0;
    for (int i = 0; i < 6; i++)
      angle += x[j].iarg[i] * fargs[i];
    // WRONG ! we must keep negative signs!
    // angle = dso::norm_angle<double, dso::AngleUnit::Radians>(angle);
    angle = std::fmod(angle, iers2010::D2PI);
    // Compute contribution from the j-th term to the polar motion coordinates
    const double sa = std::sin(angle);
    const double ca = std::cos(angle);
    dx += x[j].xs * sa + x[j].xc * ca;
    dy += x[j].ys * sa + x[j].yc * ca;
  }

  if (iband != 1) {
    const double fyears = mjd.as_fractional_years();
    // Add the secular term of the model
    dx += xrate * fyears;
    dy += yrate * fyears;
  }

  //  Finished
  return 0;
}

int iers2010::pmsdnut2(const dso::TwoPartDate &mjd, double &dx, double &dy) noexcept {
  double fargs[6];
  iers2010::utils::eop_fundarg(mjd,fargs);
  return utils::pmsdnut2(mjd,fargs,dx,dy);
}
