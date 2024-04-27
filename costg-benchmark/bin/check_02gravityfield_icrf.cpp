#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "gravity.hpp"
#include "icgemio.hpp"

constexpr const int DEGREE = 180;
constexpr const int ORDER = 180;

using namespace costg;

int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(
        stderr,
        "Usage: %s [00orbit_itrf.txt] [01earthRotation_rotaryMatrix.txt] [GEOPOTENTIAL] [02gravityfield_icrf]\n",
        argv[0]);
    return 1;
  }

  /* read orbit from input file */
  const auto orbvec = parse_orbit(argv[1]);

  /* read rotary matrix (GCRS to ITRS) from input file */
  const auto rotvec = parse_rotary(argv[2]);

  /* read accleration from input file */
  const auto accvec = parse_acceleration(argv[4]);

  /* read gravity model into a StokesCoeffs instance */
  dso::Icgem icgem(argv[3]);
  dso::StokesCoeffs stokes;
  dso::Datetime<dso::nanoseconds> t(
      dso::from_mjdepoch<dso::nanoseconds>(orbvec[0].epoch));
  if (icgem.parse_data(DEGREE, ORDER, t, stokes)) {
    fprintf(stderr, "ERROR Failed reading gravity model!\n");
    return 1;
  }

  /* checks */
  assert(stokes.max_degree() == DEGREE);
  assert(stokes.max_order() == ORDER);

  /* allocate scratch space for computations */
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> W(DEGREE + 3,
                                                                    DEGREE + 3);
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> M(DEGREE + 3,
                                                                    DEGREE + 3);

  /* spit out a title for plotting */
  printf("#title Gravity Field %s Diffs (ICRF)\n", basename(argv[3]));

  /* compare results epoch by epoch */
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 3> g;
  auto acc = accvec.begin();
  auto rot = rotvec.begin();
  for (const auto &in : orbvec) {
    /* compute acceleration for given epoch/position */
    if (dso::sh2gradient_cunningham(stokes, in.xyz, a, g, DEGREE, ORDER, -1, -1,
                                    &W, &M)) {
      fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
      return 1;
    }
    /* transform accleration vector (ITRF to GCRF) */
    const Eigen::Matrix<double,3,3> R = rot->R.transpose();
    assert(rot->epoch == in.epoch);
    a = R * a;

    /* get COSTG result */
    if (acc->epoch != in.epoch) {
      fprintf(stderr, "ERROR Faile to match epochs in input files\n");
      return 1;
    }
    printf("%d %.9f %.17e %.17e %.17e %.17e %.17e %.17e\n",
           in.epoch.imjd(), in.epoch.seconds(), acc->axyz(0), acc->axyz(1),
           acc->axyz(2), a(0), a(1), a(2));
    ++acc;
    ++rot;
  }

  return 0;
}
