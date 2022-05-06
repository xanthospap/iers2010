#include "iers2010.hpp"

// Speed of all terms given in radians per second
constexpr const double speed[] = {
    1.405190e-4, 1.454440e-4, 1.378800e-4, 1.458420e-4,
    0.729210e-4, 0.675980e-4, 0.725230e-4, 0.649590e-4,
    0.053234e-4, 0.026392e-4, 0.003982e-4,
};

// these are not used for now ...
// const double sigm2  = 1.40519e-4;
// const double sigs2  = 1.45444e-4;
// const double sign2  = 1.37880e-4;
// const double sigk2  = 1.45842e-4;
// const double sigk1  = 0.72921e-4;
// const double sigo1  = 0.67598e-4;
// const double sigp1  = 0.72523e-4;
// const double sigq1  = 0.64959e-4;
// const double sigmf  = 0.053234e-4;
// const double sigmm  = 0.026392e-4;
// const double sigssa = 0.003982e-4;

constexpr const double angfac[][4] = {
    {0.200e+01, -0.200e+01, 0.000e+00, 0.000e+00},
    {0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
    {0.200e+01, -0.300e+01, 0.100e+01, 0.000e+00},
    {0.200e+01, 0.000e+00, 0.000e+00, 0.000e+00},
    {0.100e+01, 0.000e+00, 0.000e+00, 0.250e+00},
    {0.100e+01, -0.200e+01, 0.000e+00, -0.250e+00},
    {-0.100e+01, 0.000e+00, 0.000e+00, -0.250e+00},
    {0.100e+01, -0.300e+01, 0.100e+01, -0.250e+00},
    {0.000e+00, 0.200e+01, 0.000e+00, 0.000e+00},
    {0.000e+00, 0.100e+01, -0.100e+01, 0.000e+00},
    {0.200e+01, 0.000e+00, 0.000e+00, 0.000e+00}};

int iers2010::arg2(int iyear, double day, double *angle) noexcept {
  // Constants
  constexpr const int k(11);
  constexpr const int iymin(1974);
  constexpr const double dtr(0.174532925199e-1);
  constexpr const double TWOPI(iers2010::D2PI);

  //  Validate year
  if (iyear < iymin)
    return -1;

  // Initialize day of year
  double id, fraction;
  fraction = std::modf(day, &id);

  // Compute fractional part of day in seconds
  const double fday(fraction * 86400e0);
  // Revision 07 October 2011: ICAPD modified
  const int icapd((int)id + 365 * (iyear - 1975) + ((iyear - 1973) / 4));
  const double capt((27392.500528e0 + 1.000000035e0 * (double)icapd) / 36525e0);

  // Compute mean longitude of Sun at beginning of day
  const double h0((279.69668e0 + (36000.768930485e0 + 3.03e-4 * capt) * capt) *
                  dtr);

  // Compute mean longitude of Moon at beginning of day
  const double s0(
      (((1.9e-6 * capt - .001133e0) * capt + 481267.88314137e0) * capt +
       270.434358e0) *
      dtr);

  // Compute mean longitude of lunar perigee at beginning of day
  const double p0(
      (((-1.2e-5 * capt - .010325e0) * capt + 4069.0340329577e0) * capt +
       334.329653e0) *
      dtr);

  // Compute the tidal angle arguments
  for (int i = 0; i < k; i++) {
    angle[i] = speed[i] * fday + angfac[i][0] * h0 + angfac[i][1] * s0 +
               angfac[i][2] * p0 + angfac[i][3] * TWOPI;
    angle[i] = std::fmod(angle[i], TWOPI);
    while (angle[i] < 0e0)
      angle[i] += TWOPI;
  }

  // Finished.
  return 0;
}
