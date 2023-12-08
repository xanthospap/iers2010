#include "iau.hpp"
#include "fundarg.hpp"
#ifdef USE_KAHAN_XYCIP
#include "kahan.hpp"
#endif
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

namespace {
struct Freq {
  int8_t f[5];
  int16_t nx;
  int16_t ny;
  double arg(const double *const fa) const noexcept {
    double arg = .0e0;
    for (int i = 0; i < 5; i++)
      arg += f[i] * fa[i];
    return arg;
  }
}; /* Freq */
struct XyPlntry {
  int16_t sumpl;
  int8_t f[9];
  double as;
  double ac;
  double arg(const double *const fa) const noexcept {
    if (sumpl) {
      double arg = .0e0;
      for (int i = 0; i < 9; i++)
        arg += f[i] * fa[i];
      return arg;
    }
    return 0;
  }
};
} /* unnamed namespace */

#ifdef USE_KAHAN_XYCIP
  using DSumType = dso::KahanSum;
#else
  using DSumType = double;
#endif

int dso::detail::xycip06a(const double *const fargs, double t, double &xcip,
                          double &ycip) noexcept {

  /* sum for j=4 */
  DSumType x4cip(0e0), y4cip(0e0);
  {
    /* Pointer to planetary multipliers/amplitudes for components of X and Y
     * The first nx elements are for X, the next ny components are for Y
     */
    auto it = XyPlnJ4.cbegin();
    for (const auto &s : XyLunSolJ4) {
      /* ARGUMENT computed using only luni-solar terms */
      const double arg_ls = s.arg(fargs);
      for (int i = 0; i < s.nx; i++) {
        /* add planetary contribution for X CIP */
        const double arg = arg_ls + it->arg(fargs + 5);
        const double ca = std::cos(arg);
        const double sa = std::sin(arg);
        /* add contribution from this frequency */
        x4cip += it->ac * ca + it->as * sa;
        ++it;
      }
      for (int i = 0; i < s.ny; i++) {
        /* add planetary contribution for Y CIP */
        const double arg = arg_ls + it->arg(fargs + 5);
        const double ca = std::cos(arg);
        const double sa = std::sin(arg);
        /* add contribution from this frequency (i.e. i) */
        y4cip += it->ac * ca + it->as * sa;
        ++it;
      }
    }
  }

  /* sum for j=3 */
  DSumType x3cip(0e0), y3cip(0e0);
  {
    auto it = XyPlnJ3.cbegin();
    for (const auto &s : XyLunSolJ3) {
      /* ARGUMENT computed using only luni-solar terms */
      const double arg_ls = s.arg(fargs);
      for (int i = 0; i < s.nx; i++) {
        /* add planetary contribution for X CIP */
        const double arg = arg_ls + it->arg(fargs + 5);
        const double ca = std::cos(arg);
        const double sa = std::sin(arg);
        /* add contribution from this frequency */
        x3cip += it->ac * ca + it->as * sa;
        ++it;
      }
      for (int i = 0; i < s.ny; i++) {
        /* add planetary contribution for Y CIP */
        const double arg = arg_ls + it->arg(fargs + 5);
        const double ca = std::cos(arg);
        const double sa = std::sin(arg);
        /* add contribution from this frequency (i.e. i) */
        y3cip += it->ac * ca + it->as * sa;
        ++it;
      }
    }
  }

  /* sum for j=2 */
  DSumType x2cip(0e0), y2cip(0e0);
  {
    auto it = XyPlnJ2.cbegin();
    for (const auto &s : XyLunSolJ2) {
      /* ARGUMENT computed using only luni-solar terms */
      const double arg_ls = s.arg(fargs);
      for (int i = 0; i < s.nx; i++) {
        /* add planetary contribution for X CIP */
        const double arg = arg_ls + it->arg(fargs + 5);
        const double ca = std::cos(arg);
        const double sa = std::sin(arg);
        /* add contribution from this frequency */
        x2cip += it->ac * ca + it->as * sa;
        ++it;
      }
      for (int i = 0; i < s.ny; i++) {
        /* add planetary contribution for Y CIP */
        const double arg = arg_ls + it->arg(fargs + 5);
        const double ca = std::cos(arg);
        const double sa = std::sin(arg);
        /* add contribution from this frequency (i.e. i) */
        y2cip += it->ac * ca + it->as * sa;
        ++it;
      }
    }
  }

  /* sum for j=1 */
  DSumType x1cip(0e0), y1cip(0e0);
  {
    auto it = XyPlnJ1.cbegin();
    for (const auto &s : XyLunSolJ1) {
      /* ARGUMENT computed using only luni-solar terms */
      const double arg_ls = s.arg(fargs);
      for (int i = 0; i < s.nx; i++) {
        /* add planetary contribution for X CIP */
        const double arg = arg_ls + it->arg(fargs + 5);
        const double ca = std::cos(arg);
        const double sa = std::sin(arg);
        /* add contribution from this frequency */
        x1cip += it->ac * ca + it->as * sa;
        ++it;
      }
      for (int i = 0; i < s.ny; i++) {
        /* add planetary contribution for Y CIP */
        const double arg = arg_ls + it->arg(fargs + 5);
        const double ca = std::cos(arg);
        const double sa = std::sin(arg);
        /* add contribution from this frequency (i.e. i) */
        y1cip += it->ac * ca + it->as * sa;
        ++it;
      }
    }
  }

  /* sum for j=0 */
  DSumType x0cip(0e0), y0cip(0e0);
  {
    auto it = XyPlnJ0.cbegin();
    for (const auto &s : XyLunSolJ0) {
      /* ARGUMENT computed using only luni-solar terms */
      const double arg_ls = s.arg(fargs);
      for (int i = 0; i < s.nx; i++) {
        /* add planetary contribution for X CIP */
        const double arg = arg_ls + it->arg(fargs + 5);
        const double ca = std::cos(arg);
        const double sa = std::sin(arg);
        /* add contribution from this frequency */
        x0cip += it->ac * ca + it->as * sa;
        ++it;
      }
      for (int i = 0; i < s.ny; i++) {
        /* add planetary contribution for Y CIP */
        const double arg = arg_ls + it->arg(fargs + 5);
        const double ca = std::cos(arg);
        const double sa = std::sin(arg);
        /* add contribution from this frequency (i.e. i) */
        y0cip += it->ac * ca + it->as * sa;
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
#ifdef USE_KAHAN_XYCIP
  const double x0 = (double)x0cip;
  const double x1 = (double)x1cip;
  const double x2 = (double)x2cip;
  const double x3 = (double)x3cip;
  const double x4 = (double)x4cip;
  const double y0 = (double)y0cip;
  const double y1 = (double)y1cip;
  const double y2 = (double)y2cip;
  const double y3 = (double)y3cip;
  const double y4 = (double)y4cip;
#else
  const double x0 = x0cip;
  const double x1 = x1cip;
  const double x2 = x2cip;
  const double x3 = x3cip;
  const double x4 = x4cip;
  const double y0 = y0cip;
  const double y1 = y1cip;
  const double y2 = y2cip;
  const double y3 = y3cip;
  const double y4 = y4cip;
#endif
  const auto xseries = x0 + (x1 + (x2 + (x3 + x4 * t) * t) * t) * t;
  xcip = xseries + xpoly;

  const auto yseries = y0 + (y1 + (y2 + (y3 + y4 * t) * t) * t) * t;
  ycip = yseries + ypoly;

  return 0;
}

int dso::xycip06a_new(const dso::MjdEpoch &tt, double &xcip, double &ycip,
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

  /* X-coordinate of CIP [μas] */
  /* Y-coordinate of CIP [μas] */
  dso::detail::xycip06a(fa, t, xcip, ycip);

  /* transform to [rad] and return */
  xcip = dso::sec2rad(xcip * 1e-6);
  ycip = dso::sec2rad(ycip * 1e-6);

  /* if we are given a large enough array, store the fundamental arguments */
  if (outargs)
    std::memcpy(outargs, fa, sizeof(double) * 14);

  return 0;
}
