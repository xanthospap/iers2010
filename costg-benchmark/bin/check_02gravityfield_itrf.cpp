#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "gravity.hpp"
#include "icgemio.hpp"
#include "costg_utils.hpp"

constexpr const int DEGREE = 180;
constexpr const int ORDER = 180;

using namespace costg;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(
        stderr,
        "Usage: %s [00orbit_itrf.txt] [GEOPOTENTIAL] [02gravityfield_itrf]\n",
        argv[0]);
    return 1;
  }

  /* read orbit from input file */
  const auto orbvec = parse_orbit(argv[1]);
  printf("#Note: read %d data sets from input file\n", (int)orbvec.size());

  /* read accleration from input file */
  const auto accvec = parse_acceleration(argv[3]);
  printf("#Note: read %d data sets from input file\n", (int)accvec.size());

  /* read gravity model into a StokesCoeffs instance */
  dso::Icgem icgem(argv[2]);
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
  // assert(std::abs(-0.69378736267e-09 - stokes.C(150, 150))<1e-15);

  /* allocate scratch space for computations */
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> W(DEGREE + 3,
                                                                    DEGREE + 3);
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> M(DEGREE + 3,
                                                                    DEGREE + 3);

  /* compare results epoch by epoch */
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 3> g;
  auto acc = accvec.begin();
  for (const auto &in : orbvec) {
    /* compute acceleration for given epoch/position */
    if (dso::sh2gradient_cunningham(stokes, in.xyz, a, g, DEGREE, ORDER, -1, -1, &W,
                                    &M)) {
      fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
      return 1;
    }
    /* get COSTG result */
    if (acc->epoch != in.epoch) {
      fprintf(stderr, "ERROR Faile to match epochs in input files\n");
      return 1;
    }
    printf("%d %.9f %.15e %.15e %.15e %.15e %.15e %.15e\n", in.epoch.imjd(),
           in.epoch.seconds(), acc->axyz(0), acc->axyz(1), acc->axyz(2), a(0),
           a(1), a(2));
    ++acc;
  }

  return 0;
}
