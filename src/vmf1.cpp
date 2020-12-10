#include "iers2010.hpp"
#ifdef USE_EXTERNAL_CONSTS
#include "gencon.hpp"
#endif

/// @details  This subroutine determines the Vienna Mapping Function 1 (VMF1).
///           This is the site dependent version.
///           This function is a translation/wrapper for the fortran VMF1.F
///           subroutine, found here :
///           http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  ah    Hydrostatic coefficient a (Note 1)
/// @param[in]  aw    Wet coefficient a (Note 1)
/// @param[in]  dmjd  Modified Julian Date
/// @param[in]  dlat  Ellipsoidal latitude given in radians
/// @param[in]  zd    Zenith distance in radians
/// @param[out] vmf1h Hydrostatic mapping function (Note 2)
/// @param[out] vmf1w Wet mapping function (Note 2)
/// @return           An integer, always zero
///
/// @note
///    -# The coefficients can be obtained from the website
///       http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/
///    -# The mapping functions are dimensionless scale factors
///    -# Status: Class 1 model
///
/// @version 12.01.2012
///
/// @cite iers2010
///  Boehm, J., Werl, B., and Schuh, H., (2006),
///     "Troposhere mapping functions for GPS and very long baseline
///     interferometry from European Centre for Medium-Range Weather
///     Forecasts operational analysis data," J. Geophy. Res., Vol. 111,
///     B02406, doi:10.1029/2005JB003629
///
///     Please mind that the coefficients in this paper are wrong.
///     The corrected version of the paper can be found at:
///     http://ggosatm.hg.tuwien.ac.at/DOCS/PAPERS/2006Boehm_etal_VMF1.pdf
int iers2010::vmf1(double ah, double aw, double dmjd, double dlat, double zd,
                   double &vmf1h, double &vmf1w) noexcept {
#ifdef USE_EXTERNAL_CONSTS
  constexpr double TWOPI(D2PI);
  constexpr double PI(DPI);
#else
  constexpr double PI(3.1415926535897932384626433e0);
  constexpr double TWOPI(6.283185307179586476925287e0);
#endif

  // Reference day is 28 January 1980
  // This is taken from Niell (1996) to be consistent (See References)
  const double doy{dmjd - 44239e0 + 1e0 - 28e0};

  double bh{.0029e0};
  double c0h{.062e0};

  double phh, c11h, c10h;
  if (dlat < 0e0) { // southern hemisphere
    phh = PI;
    c11h = .007e0;
    c10h = .002e0;
  } else { // northern hemisphere
    phh = 0e0;
    c11h = .005e0;
    c10h = .001e0;
  }

  const double ch{c0h + ((std::cos(doy / 365.25e0 * TWOPI + phh) + 1e0) * c11h / 2e0 +
                   c10h) *
                      (1e0 - std::cos(dlat))};

  double sine{std::sin(PI / 2e0 - zd)};
  double beta{bh / (sine + ch)};
  double gamma{ah / (sine + beta)};
  // January 10, 2012 Variable TOPCON corrected
  double topcon{(1e0 + ah / (1e0 + bh / (1e0 + ch)))};
  vmf1h = topcon / (sine + gamma);

  const double bw{0.00146e0};
  const double cw{0.04391e0};
  beta = bw / (sine + cw);
  gamma = aw / (sine + beta);
  topcon = (1e0 + aw / (1e0 + bw / (1e0 + cw)));
  vmf1w = topcon / (sine + gamma);

  // Finished
  return 0;
}
