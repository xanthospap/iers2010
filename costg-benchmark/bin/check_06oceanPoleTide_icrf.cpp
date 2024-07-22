#include "aod1b_data_stream.hpp"
#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "fundarg.hpp" /* for fundarg */
#include "geodesy/units.hpp"
#include "gravity.hpp"
#include "iau.hpp" /* for gmst */
#include "icgemio.hpp"
#include "pole_tide.hpp"

constexpr const int DEGREE = 180;
constexpr const int ORDER = 180;
constexpr const int formatD3Plot = 0;

using namespace costg;

int main(int argc, char *argv[]) {
  if (argc != 6) {
    fprintf(stderr,
            "Usage: %s [00orbit_itrf.txt] [01earthRotation_rotaryMatrix.txt] "
            "[01earthRotation_interpolatedEOP.txt] [06oceanPoleTide_icrf.txt] [desaiscopolecoef.txt]\n",
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
  const auto accvec = parse_acceleration(argv[4]);

  /* read rotary matrix (GCRS to ITRS) from input file */
  const auto rotvec = parse_rotary(argv[2]);

  /* read interpolated EOPs */
  const auto ieops = parse_eops(argv[3]);

  /* Ocean Pole tide instance */
  dso::OceanPoleTide opt(argv[5]);

  /* spit out a title for plotting */
  if (formatD3Plot) {
    printf("mjd,sec,refval,val,component\n");
  } else {
    printf("#title Pole Tide\n");
  }

  /* compare results epoch by epoch */
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 3> g;
  auto acc = accvec.begin();
  auto rot = rotvec.begin();
  auto eop = ieops.begin();
  for (const auto &in : orbvec) {
    /* GPSTime */
    const auto tt = in.epoch;

    /* (xp,yp) rad to arcsec */
    const double xp = dso::rad2sec(eop->xp);
    const double yp = dso::rad2sec(eop->yp);
    assert(eop->epoch == in.epoch);

    /* compute potential corrections from Pole Tide */
    if (opt.stokes_coeffs(tt, xp, yp, DEGREE, ORDER)) {
      fprintf(stderr, "ERROR Failed computing Stokes COefficients\n");
      return 1;
    }

    /* compute acceleration for given epoch/position */
    if (dso::sh2gradient_cunningham(opt.stokes_coeffs(), in.xyz, a, g, DEGREE,
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
      printf("%d,%.9f,%.17e,%.17e,X\n", in.epoch.imjd(), in.epoch.seconds(),
             acc->axyz(0), a(0));
      printf("%d,%.9f,%.17e,%.17e,Y\n", in.epoch.imjd(), in.epoch.seconds(),
             acc->axyz(1), a(1));
      printf("%d,%.9f,%.17e,%.17e,Z\n", in.epoch.imjd(), in.epoch.seconds(),
             acc->axyz(2), a(2));
    } else {
      printf("%d %.9f %.17e %.17e %.17e %.17e %.17e %.17e\n", in.epoch.imjd(),
             in.epoch.seconds(), acc->axyz(0), acc->axyz(1), acc->axyz(2), a(0),
             a(1), a(2));
    }

    ++acc;
    ++rot;
    ++eop;
  }

  return 0;
}
