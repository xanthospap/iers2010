#include "eop.hpp"
#include "fundarg.hpp"
#include "iau.hpp"

/* Test case provided at the header of PMSDNUT2
 */

int main() {
  
  dso::MjdEpoch t(54335); /* August 23, 2007 */

  /* compute fundamental arguments */
  double fargs[6];
  dso::fundarg(t, fargs);

  /* compute gmst */
  const double dut1 = -0.1631146e0; /* UT1 - UTC */
  const auto ut1 = t.tt2ut1(dut1);
  const double gmst = dso::gmst(t, ut1);

  /* compute libration on pole */
  double xp,yp;
  dso::xy_libration(fargs, gmst, xp, yp);

  /* reference results */
  const double pm1 = 24.83144238273364834e0; /* microarcseconds */
  const double pm2 = -14.09240692041837661e0; /* microarcseconds */

  printf("delta xp %+.6f diff %+.6f [GMST06]\n", pm1, pm1-xp);
  printf("delta yp %+.6f diff %+.6f\n", pm2, pm2-yp);

  const double era = dso::era00(ut1);
  dso::xy_libration(fargs, era, xp, yp);
  printf("delta xp %+.6f diff %+.6f [ERA]\n", pm1, pm1-xp);
  printf("delta yp %+.6f diff %+.6f\n", pm2, pm2-yp);
  
  const double g82 = dso::gmst82(ut1);
  dso::xy_libration(fargs, g82, xp, yp);
  printf("delta xp %+.6f diff %+.6f [GMST82]\n", pm1, pm1-xp);
  printf("delta yp %+.6f diff %+.6f\n", pm2, pm2-yp);

  double tc = t.jcenturies_sinceJ2000();
  double gmst_old =
      std::fmod(67310.54841e0 + tc * ((8640184.812866e0 + 3155760000e0) +
                                     tc * (0.093104e0 + tc * (-0.0000062))),
                86400e0);
  gmst_old = gmst_old / (86400e0/dso::D2PI);
  dso::xy_libration(fargs, gmst_old, xp, yp);
  printf("delta xp %+.6f diff %+.6f [IERS10]\n", pm1, pm1-xp);
  printf("delta yp %+.6f diff %+.6f\n", pm2, pm2-yp);
  
  return 0;
}
