#include "fundarg.hpp"
#include "iau.hpp"
#include "foo.hpp"
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

int dso::detail::xycip06a(const double *const fargs, double t, double &xcip,
                          double &ycip) noexcept {
  /* powers of t */
  double tpow[5] = {1e0};
  for (int i = 1; i < 5; i++)
    tpow[i] = t * tpow[i - 1];

  /* sum for j=0,..4 */
  double xser(0e0), yser(0e0);
  {
    auto it = XyPln.cbegin();
    // const XyPlntry *__restrict__ it = XyPln.data();
    for (const auto &s : XyLunSol) {
      /* ARGUMENT computed using only luni-solar terms */
      const double arg_ls = s.arg(fargs);
      const double ca = std::cos(arg_ls);
      const double sa = std::sin(arg_ls);
      for (int i = 0; i < s.nx; i++) {
        /* add planetary contribution for X CIP */
        const int plterms = it->sumpl;
        const double arg = plterms ? (arg_ls + it->arg(fargs + 5)) : arg_ls;
        const double cax = plterms ? std::cos(arg) : ca;
        const double sax = plterms ? std::sin(arg) : sa;
        /* add contribution from this frequency */
        xser += (it->ac * cax + it->as * sax) * tpow[it->pow];
        ++it;
      }
      for (int i = 0; i < s.ny; i++) {
        /* add planetary contribution for Y CIP */
        const int plterms = it->sumpl;
        const double arg = plterms ? (arg_ls + it->arg(fargs + 5)) : arg_ls;
        const double cay = plterms ? std::cos(arg) : ca;
        const double say = plterms ? std::sin(arg) : sa;
        /* add contribution from this frequency (i.e. i) */
        yser += (it->ac * cay + it->as * say) * tpow[it->pow];
        ++it;
      }
    }
  }

  /* Add polynomial terms in [microarcseconds] */
  constexpr const double XPolyCoeffs[] = {
      -16617e0, 2004191898e0, -429782.9e0, -198618.34e0, 7.578e0, 5.9285e0};
  const double xpoly =
      XPolyCoeffs[0] +
      (XPolyCoeffs[1] +
       (XPolyCoeffs[2] +
        (XPolyCoeffs[3] + (XPolyCoeffs[4] + XPolyCoeffs[5] * t) * t) * t) *
           t) *
          t;

  /* Add polynomial terms in [microarcseconds] */
  constexpr const double YPolyCoeffs[] = {-6951e0,   -25896e0,   -22407274.7e0,
                                          1900.59e0, 1112.526e0, 0.1358e0};
  const double ypoly =
      YPolyCoeffs[0] +
      (YPolyCoeffs[1] +
       (YPolyCoeffs[2] +
        (YPolyCoeffs[3] + (YPolyCoeffs[4] + YPolyCoeffs[5] * t) * t) * t) *
           t) *
          t;

  /* accumulate and return in [microarcseconds] */
  xcip = xser + xpoly;
  ycip = yser + ypoly;

  return 0;
}

int dso::xycip06a(const dso::MjdEpoch &tt, double &xcip, double &ycip,
                  double *outargs) noexcept {
  /* Interval between fundamental date J2000.0 and given date. */
  const double t = tt.jcenturies_sinceJ2000();

  /* Fundamental arguments (IERS 2003) */
  const double fa[14] = {
      dso::iers2010::fal03(t),  dso::iers2010::falp03(t),
      dso::iers2010::faf03(t),  dso::iers2010::fad03(t),
      dso::iers2010::faom03(t), dso::iers2010::fame03(t),
      dso::iers2010::fave03(t), dso::iers2010::fae03(t),
      dso::iers2010::fama03(t), dso::iers2010::faju03(t),
      dso::iers2010::fasa03(t), dso::iers2010::faur03(t),
      dso::iers2010::fane03(t), dso::iers2010::fapa03(t),
  };

  /* XY -coordinates of CIP [μas] */
  dso::detail::xycip06a(fa, t, xcip, ycip);

  /* transform to [rad] and return */
  xcip = dso::sec2rad(xcip * 1e-6);
  ycip = dso::sec2rad(ycip * 1e-6);

  /* if we are given a large enough array, store the fundamental arguments */
  if (outargs)
    std::memcpy(outargs, fa, sizeof(double) * 14);

  return 0;
}
