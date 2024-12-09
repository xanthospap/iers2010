#include "aod1b_data_stream.hpp"
#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "gravity.hpp"
#include "icgemio.hpp"

constexpr const int DEGREE = 180;
constexpr const int ORDER = 180;
constexpr const int formatD3Plot = 0;

using namespace costg;

int main(int argc, char *argv[]) {
  if (argc != 4 && argc != 6) {
    fprintf(stderr,
            "Usage: Usage: %s [00orbit_itrf.txt] "
            "[AOD1B_2008-07-03_X_06.asc] [08aod1b_RL06_icrf.txt] "
            "[01earthRotation_rotaryMatrix.txt] [AOD1B_DATA_DIR]\n",
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

  /* read acceleration from input file */
  const auto accvec = parse_acceleration(argv[3]);

  std::vector<BmRotaryMatrix> rotvec;
  /* read rotary matrix (GCRS to ITRS) from input file */
  rotvec = parse_rotary(argv[4]);

  /* we initialize the instance using a starting (AOD1B) file and a directory
   * where subsequent files are placed and can be used.
   */
  dso::Aod1bDataStream<dso::AOD1BCoefficientType::GLO> aodin(argv[2], argv[5]);
  aodin.initialize();

  dso::StokesCoeffs stokes(DEGREE, ORDER, aodin.stream().GM(),
                           aodin.stream().Re());

  /* spit out a title for plotting */
  if (formatD3Plot) {
    printf("mjd,sec,refval,val,component\n");
  } else {
    printf("#title De-aliasing (data: %s)\n", basename(argv[2]));
  }

  /* compare results epoch by epoch */
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 3> g;
  auto acc = accvec.begin();
  auto rot = rotvec.begin();
  for (const auto &in : orbvec) {
    /* GPSTime */
    const auto t = in.epoch;
    dso::datetime<dso::nanoseconds> dt(dso::from_mjdepoch<dso::nanoseconds>(t));

    /* get Stokes coefficients for this epoch from the AOD1B file */
    if (aodin.coefficients_at(dt.gps2tt(), stokes)) {
      fprintf(stderr, "Failed interpolating coefficients\n");
      return 1;
    }

    /* for the test, degree one coefficients are not taken into account */
    stokes.C(0, 0);
    stokes.C(1, 0) = stokes.C(1, 1) = 0e0;
    stokes.S(1, 1) = 0e0;

    /* compute acceleration for given epoch/position */
    if (dso::sh2gradient_cunningham(stokes, in.xyz, a, g, DEGREE, ORDER, -1, -1,
                                    &W, &M)) {
      fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
      return 1;
    }

    /* if needed, transform acceleration from ITRF to GCRF */
    const Eigen::Matrix<double, 3, 3> R = rot->R.transpose();
    assert(rot->epoch == in.epoch);
    a = R * a;

    /* get COSTG result */
    if (acc->epoch != in.epoch) {
      fprintf(stderr, "ERROR Faile to match epochs in input files\n");
      return 1;
    }

    if (formatD3Plot) {
      printf("%d,%.9f,%.17e,%.17e,X\n", in.epoch.imjd(),
             in.epoch.seconds().seconds(), acc->axyz(0), a(0));
      printf("%d,%.9f,%.17e,%.17e,Y\n", in.epoch.imjd(),
             in.epoch.seconds().seconds(), acc->axyz(1), a(1));
      printf("%d,%.9f,%.17e,%.17e,Z\n", in.epoch.imjd(),
             in.epoch.seconds().seconds(), acc->axyz(2), a(2));
    } else {
      printf("%d %.9f %.17e %.17e %.17e %.17e %.17e %.17e\n", in.epoch.imjd(),
             in.epoch.seconds().seconds(), acc->axyz(0), acc->axyz(1),
             acc->axyz(2), a(0), a(1), a(2));
    }

    ++acc;
    ++rot;
  }

  return 0;
}
