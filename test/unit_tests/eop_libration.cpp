#include "eop.hpp"
#include "fundarg.hpp"
#include "iau.hpp"
#include <charconv>
#include <fstream>

constexpr const double MAX_MICROARCSEC = 1e0;
[[maybe_unused]] constexpr const double MAX_MICROSEC = 1e0;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "USAGE: %s [DATA]\n", argv[0]);
    return 1;
  }

  std::ifstream fin(argv[1]);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed opening data file %s\n", argv[1]);
    fprintf(stderr,
            "Note that reference results for this program can be produced via "
            "the test/fortran/eop-libration.out program\n");
    return 2;
  }

  char line[512];
  double d[5];
  int error = 0;
  while (fin.getline(line, 512) && (!error)) {

    int sz = std::strlen(line);
    const char *s = line;
    for (int i = 0; i < 5; i++) {
      while (*s && *s == ' ')
        ++s;
      auto t = std::from_chars(s, line + sz, d[i]);
      if (t.ec != std::errc{})
        ++error;
      s = t.ptr;
    }

    /* Setup date in mjd */
    int mjd = (int)d[0];
    double fdaysec = (d[0] - mjd) * 86400e0;
    dso::MjdEpoch t(mjd, dso::FractionalSeconds(fdaysec));

    /* compute fundamental arguments */
    double fargs[6];
    dso::fundarg(t, fargs);

    /* compute approximate gmst */
    const double gmst = dso::gmst82(t);

    /* compute libration variations on xp, yp */
    double xp, yp, dut1, dlod;
    dso::deop_libration(fargs, gmst, xp, yp, dut1, dlod);

    // printf("[MINE] %20.5f %+.6f %+.6f %+.6f %+.6f\n", fdaysec, xp, yp, dut1,
    // dlod); printf("[IERS] %20.5f %+.6f %+.6f %+.6f %+.6f\n", fdaysec, d[1],
    // d[2],d[3],d[4]);

    assert(std::abs(d[1] - xp) < MAX_MICROARCSEC);
    assert(std::abs(d[2] - yp) < MAX_MICROARCSEC);
    assert(std::abs(d[3] - dut1) < MAX_MICROSEC);
    assert(std::abs(d[4] - dlod) < MAX_MICROSEC);
  }

  return 0;
}
