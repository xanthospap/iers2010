#include "iersconst.hpp"
#include "pole_tide.hpp"
#include "geodesy/units.hpp"

int dso::oceanPoleTide::stokes_coeffs(const dso::MjdEpoch &t, double xp,
                                      double yp, int max_degree, int max_order,
                                      dso::StokesCoeffs &cs) noexcept {

  if ((max_degree < max_order) ||
      (max_degree > dso::oceanPoleTide::MAX_DEGREE_DESAI_2002) ||
      (max_order > dso::oceanPoleTide::MAX_ORDER_DESAI_2002)) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for parsing Ocean Pole Tide %s "
            "(traceback: %s)\n",
            fn, __func__);
    return 1;
  }

  constexpr const double rhow = 1025e0; /* density of sea water in [kgm^−3] */
  constexpr const double g2_real = 0.6870e0;
  constexpr const double g2_imag = 0.0036e0;
  constexpr const double fac = ((iers2010::OmegaEarth * iers2010::OmegaEarth) *
                                std::pow(iers2010::Re, 4) / iers2010::GM) *
                               (2e0 * D2PI * iers2010::Ge * rhow) /
                               iers2010::ge;
  const auto m12 = dso::pole_tide_details::mcoeffs(t, xp, yp);
  const double m1 = m12.m1; /* m1 in [arcsec] */
  const double m2 = m12.m2; /* m2 in [arcsec] */
  const double freal = dso::sec2rad(m1 * g2_real + m2 * g2_imag);
  const double fimag = dso::sec2rad(m2 * g2_real - m1 * g2_imag);
  double knfac[MAX_DEGREE_DESAI_2002];

  /* order: m = 0; compute kn factors */
  int m = 0;
  for (int l = m; l <= max_degree; l++) {
    knfac[l] = (1e0 + kn[l]) / (2 * l + 1);
    cs.C(l, m) = fac * knfac[l] * (A_real(l, m) * freal + A_imag * fimag);
    cs.S(l, m) = fac * knfac[l] * (B_real(l, m) * freal + B_imag * fimag);
  }

  /* all toher orders (m) except m=0 */
  for (m = 1; m <= max_order; m++) {
    for (int l = m; l <= max_degree; l++) {
      cs.C(l, m) = fac * knfac[l] * (A_real(l, m) * freal + A_imag * fimag);
      cs.S(l, m) = fac * knfac[l] * (B_real(l, m) * freal + B_imag * fimag);
    }
  }

  return 0;
}
