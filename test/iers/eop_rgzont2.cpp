
#include "eop.hpp"
#include "fundarg.hpp"
#include "iau.hpp"
#include <charconv>
#include <fstream>
#ifdef NDEBUG
#undef NDEBUG
#endif

constexpr const double MAX_MICROSEC = 1e4;
constexpr const double MAX_MICROSECDAY = 1e+3;
constexpr const double MAX_RADSEC = 1e+0;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "USAGE: %s [DATA]\n", argv[0]);
    fprintf(stderr, "Note that reference results for this program can be "
                    "produced via the test/fortran/TEST_RG_ZONT2 program\n");
    return 1;
  }

  std::ifstream fin(argv[1]);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed opening data file %s\n", argv[1]);
    return 2;
  }

  /* read input values (epoch, dut1, dlod, domega) in sec, sec/day and rad/sec
   */
  char line[512];
  double d[4];
  int error = 0;
  while (fin.getline(line, 512) && (!error)) {

    int sz = std::strlen(line);
    const char *s = line;
    for (int i = 0; i < 4; i++) {
      while (*s && *s == ' ')
        ++s;
      auto t = std::from_chars(s, line + sz, d[i]);
      if (t.ec != std::errc{})
        ++error;
      s = t.ptr;
    }
    if (error) {
      fprintf(stderr, "ERROR. Failed parsing input data!\n");
      assert(1 == 2);
    }

    /* Setup date in mjd */
    int mjd = (int)d[0];
    double fdaysec = (d[0] - mjd) * 86400e0;
    dso::MjdEpoch t(mjd, dso::FractionalSeconds(fdaysec));

    /* compute fundamental arguments */
    double fargs[6];
    dso::fundarg(t, fargs);

    /* compute ocean tide variations on EOPs */
    double dut1, dlod, domega;
    dso::deop_zonal_tide(fargs, dut1, dlod, domega);

    /* scale */
    dut1 *= 1e+6;
    dlod *= 1e+6;
    domega *= 1e+14;

    printf("%20.5f (%.3f %.3f) (%.3f %.3f) %+.6f %+.6f %+.15f\n", fdaysec, d[1],
           dut1, d[2], dlod, std::abs(d[1] - dut1), std::abs(d[2] - dlod),
           std::abs(d[3] - domega));

    /* note that results are in microarcsec/microsec */
    if (!(std::abs(d[1] - dut1) * 1e6 < MAX_MICROSEC))
      std::exit(EXIT_FAILURE);
    if (!(std::abs(d[2] - dlod) * 1e6 < MAX_MICROSECDAY))
      std::exit(EXIT_FAILURE);
    if (!(std::abs(d[3] - domega) < MAX_RADSEC))
      std::exit(EXIT_FAILURE);
  }

  return 0;
}
