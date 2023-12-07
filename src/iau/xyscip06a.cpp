#include "iau.hpp"
#include "fundarg.hpp"
#include "geodesy/units.hpp"
#include <cstring>
#include <thread>
#include <utility>

namespace {
/* Compute S <- s + xcip * ycip / 2e0 in [rad]; result (S) is stored in sxy2 */
inline void t_s06(const double *const fargs, double t,
                    double &sxy2) noexcept {
  sxy2 = dso::detail::s06(fargs, t, 0e0, 0e0);
}

/* X-coordinate of CIP [rad]; result is stored in xcip */
inline void t_xcip06a(const double *const fargs, double t,
                      double &xcip) noexcept {
  xcip = dso::sec2rad(dso::detail::xcip06a(fargs, t) * 1e-6);
}

/* Y-coordinate of CIP [rad]; result is stored in ycip */
inline void t_ycip06a(const double *const fargs, double t,
                      double &ycip) noexcept {
  ycip = dso::sec2rad(dso::detail::xcip06a(fargs, t) * 1e-6);
}
}/* unnamed namespace */

int dso::xys06a(const dso::MjdEpoch &tt, double &xcip,
                  double &ycip, double &s, double *outargs) noexcept {
  /* Interval between fundamental date J2000.0 and given date. */
  const double t = tt.jcenturies_sinceJ2000();

  /* Fundamental arguments (IERS 2003) */
  const double fa[14] = {
      dso::iers2010::fal03(t),
      dso::iers2010::falp03(t),
      dso::iers2010::faf03(t),
      dso::iers2010::fad03(t),
      dso::iers2010::faom03(t),
      dso::iers2010::fame03(t),
      dso::iers2010::fave03(t),
      dso::iers2010::fae03(t),
      dso::iers2010::fama03(t),
      dso::iers2010::faju03(t),
      dso::iers2010::fasa03(t),
      dso::iers2010::faur03(t),
      dso::iers2010::fane03(t),
      dso::iers2010::fapa03(t),
  };

  /* X-coordinate of CIP [rad] */
  std::thread t1(t_xcip06a, fa, t, std::ref(xcip));

  /* Y-coordinate of CIP [rad] */
  std::thread t2(t_ycip06a, fa, t, std::ref(ycip));

  /* This is **NOT** the CIO, but S <- s + xcip * ycip / 2e0 */
  std::thread t3(t_s06, fa, t, std::ref(s)); 

  /* wait untill all threads are over */
  t1.join();
  t2.join();
  t3.join();
  
  /* compute CIO locator: s = S - xcip * ycip / 2e0 */
  s = s - xcip*ycip/2e0;

  /* if we are given a large enough array, store the fundamental arguments */
  if (outargs) 
    std::memcpy(outargs, fa, sizeof(double) * 14);

  return 0;
}
