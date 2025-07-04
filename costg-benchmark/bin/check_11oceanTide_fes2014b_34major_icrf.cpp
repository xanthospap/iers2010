#include "ocean_tide.hpp"
#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "eop.hpp"
#include "fundarg.hpp"
#include "gravity.hpp"
#include "icgemio.hpp"

constexpr const int DEGREE = 120;
constexpr const int ORDER = 120;
constexpr const int formatD3Plot = 0;

using namespace costg;

int main(int argc, char *argv[]) {
  if (argc != 7) {
    fprintf(stderr,
            "Usage: %s [00orbit_itrf.txt] [11oceanTide_fes2014b_34major_icrf.txt] "
            "[01earthRotation_rotaryMatrix.txt] [eopc04.1962-now] "
            "[FES2014b_OCN_001fileList.txt] [FES2014b gfc dir]\n",
            argv[0]);
    return 1;
  }

  /* allocate scratch space for computations */
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> W(DEGREE + 3,
                                                                    DEGREE + 3);
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> M(DEGREE + 3,
                                                                    DEGREE + 3);

  /* read orbit from input file */
  const auto orbvec = parse_orbit(argv[1]);

  /* read accleration from input file */
  const auto accvec = parse_acceleration(argv[2]);

  /* read rotary matrix (GCRS to ITRS) from input file */
  std::vector<BmRotaryMatrix> rotvec = parse_rotary(argv[3]);

  /* read EOPS (we will need DUT1) */
  dso::EopSeries eop;
  {
    auto imjd = orbvec[0].epoch.imjd();
    const auto t1 = dso::MjdEpoch(imjd - 2);
    const auto t2 = dso::MjdEpoch(imjd + 3);
    /* parse EOP values to the EopSeries instance for any t >= t1 and t < t2 */
    if (dso::parse_iers_C04(argv[4], t1, t2, eop)) {
      fprintf(stderr, "ERROR Failed parsing eop file\n");
      return 1;
    }
  }

  /* create an OceanTide instance using the fes2014b gfc files */
  dso::OceanTide fes14b =
      dso::ocean_tide_from_gfc(argv[5], argv[6], DEGREE, ORDER);

  /* spit out a title for plotting */
  if (formatD3Plot) {
    printf("mjd,sec,refval,val,component\n");
  } else {
    printf("#title Ocean Tidal Loading - %s Major\n", fes14b.atlas().name());
  }

  /* compare results epoch by epoch */
  double fargs[5];
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 3> g;
  auto acc = accvec.begin();
  auto rot = rotvec.begin();
  for (const auto &in : orbvec) {
    /* GPSTime */
    const auto t = in.epoch;

    /* compute gmst using an approximate value for UT1 (linear interpolation) */
    double dut1_approx;
    eop.approx_dut1(t.gps2tai().tai2tt(), dut1_approx);

    /* compute fundamental (Delaunay) arguments fot t */
    dso::fundarg(t.gps2tai().tai2tt(), fargs);

    /* compute Stokes coeffs (for atm. tides) */
    fes14b.stokes_coeffs(t.gps2tai().tai2tt(), t.tt2ut1(dut1_approx), fargs);

    /* for the test, degree one coefficients are not taken into account */
    fes14b.stokes_coeffs().C(0, 0) = fes14b.stokes_coeffs().C(1, 0) = fes14b.stokes_coeffs().C(1, 1) = 0e0;
    fes14b.stokes_coeffs().S(1, 1) = 0e0;

    /* compute acceleration for given epoch/position (ITRF) */
    if (dso::sh2gradient_cunningham(fes14b.stokes_coeffs(), in.xyz, a, g, DEGREE,
                                    ORDER, -1, -1, &W, &M)) {
      fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
      return 1;
    }

    /* transform acceleration from ITRF to GCRF */
    const Eigen::Matrix<double, 3, 3> R = rot->R.transpose();
    assert(rot->epoch == in.epoch);
    a = R * a;

    /* get COSTG result */
    if (acc->epoch != in.epoch) {
      fprintf(stderr, "ERROR Failed to match epochs in input files\n");
      return 1;
    }

    if (formatD3Plot) {
      printf("%d,%.9f,%.17e,%.17e,X\n", in.epoch.imjd(), in.epoch.seconds().seconds(),
             acc->axyz(0), a(0));
      printf("%d,%.9f,%.17e,%.17e,Y\n", in.epoch.imjd(), in.epoch.seconds().seconds(),
             acc->axyz(1), a(1));
      printf("%d,%.9f,%.17e,%.17e,Z\n", in.epoch.imjd(), in.epoch.seconds().seconds(),
             acc->axyz(2), a(2));
    } else {
      printf("%d %.9f %.17e %.17e %.17e %.17e %.17e %.17e\n", in.epoch.imjd(),
             in.epoch.seconds().seconds(), acc->axyz(0), acc->axyz(1), acc->axyz(2), a(0),
             a(1), a(2));
    }

    /* Debugging */
    // printf("Cnm(5,5)sin=%.12e Snm(5,5)sin=%.12e\n", fes14b.atlas().waves()[0].stokes_sin().C(5,5), fes14b.atlas().waves()[0].stokes_sin().S(5,5));

    ++acc;
    ++rot;
  }

  return 0;
}
