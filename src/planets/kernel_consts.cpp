#include "planets.hpp"

// see https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html
int dso::sun_moon_gm(double &GMSun, double &GMMoon, int use_si,
                     const char *pck_kernel) noexcept {
  if (pck_kernel)
    dso::cspice::load_if_unloaded(pck_kernel);

  int n;
  /* [km^3/ sec^2] */
  bodvrd_c("SUN", "GM", 1, &n, &GMSun);
  bodvrd_c("MOON", "GM", 1, &n, &GMMoon);

  if (use_si) {
    GMSun *= 1e9;
    GMMoon *= 1e9;
  }

  return 0;
}
