#include "eop.hpp"
#include "fundarg.hpp"
#include "iau.hpp"
#include <fstream>
#include <charconv>

constexpr const double MAX_MICROSEC = 1e-1;
constexpr const double MAX_DOMEGA = 1e-16; /* rad/sec */

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "USAGE: %s [DATA]\n", argv[0]);
    return 1;
  }

  std::ifstream fin(argv[1]);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed opening data file %s\n", argv[1]);
    fprintf(stderr, "Note that reference results for this program can be produced via the test/fortran/test-rgzont2.out program\n");
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

    /* compute zonal tides effect on Earth's rotation */
    double dut1,dlod,domega;
    dso::deop_zonal_tide(fargs,dut1,dlod,domega);

    /* scale to compare results */
    d[1] *= 1e6; /* microseconds */
    d[2] *= 1e6; /* microseconds */
    domega *= 1e-14; /* rad/sec */

    // printf("%20.5f %+.6f %+.6f %+.16f \n", fdaysec,dut1,dlod,domega);
    // printf("                     %+.6f %+.6f %+.16f\n", d[1],d[2],d[3]);

    assert(std::abs(d[1] - dut1) < MAX_MICROSEC);
    assert(std::abs(d[2] - dlod) < MAX_MICROSEC);
    assert(std::abs(d[3] - domega) < MAX_DOMEGA);
  }

  return 0;
}
