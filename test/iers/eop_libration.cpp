#include "eop.hpp"
#include "fundarg.hpp"
#include "iau.hpp"
#include <charconv>
#include <fstream>
#ifdef NDEBUG
#undef NDEBUG
#include <cassert>
#endif

constexpr const double MAX_MICROSEC = 1e-1;
constexpr const double MAX_MICROSECDAY = 1e+0;
constexpr const double MAX_MICROARCSEC = 1e+0;

int main(int argc, char* argv[])
{
	if (argc != 2) {
		fprintf(stderr, "USAGE: %s [DATA]\n", argv[0]);
		fprintf(stderr, "Note that reference results for this program can be "
				"produced via the test/fortran/TEST_UTLIBR program\n");
		return 1;
	}

	std::ifstream fin(argv[1]);
	if (!fin.is_open()) {
		fprintf(stderr, "ERROR Failed opening data file %s\n", argv[1]);
		return 2;
	}

	/* read input values (epoch, dut1, dlod, dx, dy) in micorsec, microsec/day
	 * and microarcsec
	 */
	char line[512];
	double d[5];
	int error = 0;
	while (fin.getline(line, 512) && (!error)) {

		int sz = std::strlen(line);
		const char* s = line;
		for (int i = 0; i < 5; i++) {
			while (*s && *s == ' ')
				++s;
			auto t = std::from_chars(s, line + sz, d[i]);
			if (t.ec != std::errc {})
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

		const double gmst = dso::gmst82(t.tt2ut1(0e0));

		double dxp, dyp, dut1, dlod;
		dso::deop_libration(fargs, gmst, dxp, dyp, dut1, dlod);

		//printf("%.9f \n", dut1);
		//printf("%.9f \n", d[1]);

		/* note that results are in microarcsec/microsec */
		if (!(std::abs(d[1] - dut1) < MAX_MICROSEC))
			std::exit(EXIT_FAILURE);
		if (!(std::abs(d[2] - dlod) < MAX_MICROSECDAY))
			std::exit(EXIT_FAILURE);
		if (!(std::abs(d[3] - dxp) < MAX_MICROARCSEC))
			std::exit(EXIT_FAILURE);
		if (!(std::abs(d[4] - dyp) < MAX_MICROARCSEC))
			std::exit(EXIT_FAILURE);
	}

	return 0;
}
