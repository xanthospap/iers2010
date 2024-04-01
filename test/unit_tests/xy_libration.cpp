#include "eop.hpp"
#include "fundarg.hpp"
#include "iau.hpp"
#include <fstream>
#include <charconv>

/* Input data for the tests are produced by the source file 
 * fortran/pmsdnut2.cpp
 */

constexpr const double MAX_MICROARCSEC = 1e-6;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "USAGE: %s [PMSDNUT2 DATA]\n", argv[0]);
    return 1;
  }

  std::ifstream fin(argv[1]);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed opening PMSDNUT2 data file %s\n", argv[1]);
    return 2;
  }

  char line[512];
  double d[3];
  int error = 0;
  while (fin.getline(line, 512) && (!error)) {

    int sz = std::strlen(line);
    const char *s = line;
    for (int i=0; i<3; i++) {
      while (*s && *s == ' ') ++s;
      auto t = std::from_chars(s, line+sz, d[i]);
      if (t.ec != std::errc{}) ++error;
      s = t.ptr;
    }

    /* Setup date in mjd */
    int mjd = (int)d[0];
    double fdaysec = (d[0]-mjd)*86400e0;
    dso::MjdEpoch t(mjd, dso::FractionalSeconds(fdaysec)); 

    /* compute fundamental arguments */
    double fargs[6];
    dso::fundarg(t, fargs);

    /* compute approximate gmst */
    const double gmst = dso::gmst82(t);

    /* compute libration on pole */
    double xp, yp;
    dso::xypole_libration(fargs, gmst, xp, yp);

    assert(std::abs(d[1] - xp) < MAX_MICROARCSEC);
    assert(std::abs(d[2] - yp) < MAX_MICROARCSEC);
    //printf("delta xp %+.6f diff %+.6f\n", d[1], d[1] - xp);
    //printf("delta yp %+.6f diff %+.6f\n", d[2], d[2] - yp);
  }

  return 0;
}
