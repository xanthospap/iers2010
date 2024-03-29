#include "eop.hpp"
#include "fundarg.hpp"

/* Test case provided at the header of PMSDNUT2
 */

int main() {
  
  dso::MjdEpoch t(54335); /* August 23, 2007 */

  /* compute fundamental arguments */
  double fargs[6];
  dso::fundarg(t, fargs);

  /* compute gmst */

  /* compute libration on pole */
  double xp,yp;
  dso::xy_libration(const double *const fargs, double gmst, xp,
                 yp);
}
