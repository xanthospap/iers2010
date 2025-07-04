#include "eop.hpp"
#include "fundarg.hpp"
#include "iau.hpp"
#include <fstream>
#include <charconv>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "sofa.h"
#include <cassert>

constexpr const double MAX_MICROARCSEC = 2e-0;
constexpr const double MAX_MICROSEC = 1e-0;

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    fprintf(stderr, "USAGE: %s [DATA FILE]\n", argv[0]);
    fprintf(stderr, "Note that reference results for this program can be produced via the test/fortran/TEST_ORTHO_EOP program\n");
    return 1;
  }

  std::ifstream fin(argv[1]);
  if (!fin.is_open())
  {
    fprintf(stderr, "ERROR Failed opening data file %s\n", argv[1]);
    return 2;
  }

  char line[512];
  double d[4];
  int error = 0;
  while (fin.getline(line, 512) && (!error))
  {

    int sz = std::strlen(line);
    const char *s = line;
    for (int i = 0; i < 4; i++)
    {
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

    /* compute ocean tide variations on EOPs */
    double xp, yp, dut1, dlod;
    dso::deop_ocean_tide(fargs, gmst, xp, yp, dut1, dlod);

    assert(std::abs(d[1] - xp) < MAX_MICROARCSEC);
    assert(std::abs(d[2] - yp) < MAX_MICROARCSEC);
    assert(std::abs(d[3] - dut1) < MAX_MICROSEC);
  }

  return 0;
}
