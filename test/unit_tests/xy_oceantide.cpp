#include "eop.hpp"
#include "fundarg.hpp"
#include "iau.hpp"
#include <fstream>
#include <charconv>

/* Input data for the tests are produced by the source file 
 * fortran/ortho_eop.cpp
 */

constexpr const double MAX_MICROARCSEC = 5e-1;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "USAGE: %s [ORTHO_EOP DATA]\n", argv[0]);
    return 1;
  }

  std::ifstream fin(argv[1]);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed opening ORTHO_EOP data file %s\n", argv[1]);
    return 2;
  }

  char line[512];
  double d[4];
  int error = 0;
  while (fin.getline(line, 512) && (!error)) {

    int sz = std::strlen(line);
    const char *s = line;
    for (int i=0; i<4; i++) {
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

    /* compute oceantide effect on pole */
    double xp, yp;
    dso::xypole_oceantide(fargs, gmst, xp, yp);

    /* compute oceantide effect on UT1 and LOD */
    double dut1, dlod;
    dso::utlod_oceantide(fargs, gmst, dut1, dlod);


    printf("delta xp   %+.6f diff %+.6f\n", d[1], d[1] - xp);
    printf("delta yp   %+.6f diff %+.6f\n", d[2], d[2] - yp);
    printf("delta dut1 %+.6f diff %+.6f\n", d[3], d[3] - dut1);
    assert(std::abs(d[1] - xp) < MAX_MICROARCSEC);
    assert(std::abs(d[2] - yp) < MAX_MICROARCSEC);
    assert(std::abs(d[3] - dut1) < MAX_MICROARCSEC);
  }

  return 0;
}
