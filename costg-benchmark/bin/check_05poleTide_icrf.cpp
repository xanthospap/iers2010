#include "aod1b_data_stream.hpp"
#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "gravity.hpp"
#include "icgemio.hpp"
#include "pole_tide.hpp"
#include "fundarg.hpp" /* for fundarg */
#include "iau.hpp"     /* for gmst */
#include "geodesy/units.hpp"

constexpr const int DEGREE = 2;
constexpr const int ORDER = 2;
constexpr const int formatD3Plot = 1;

using namespace costg;

int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(stderr,
            "Usage: %s [00orbit_itrf.txt] [01earthRotation_rotaryMatrix.txt] [eopc04.1962-now] [05poleTide_icrf.txt]\n",
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

  std::vector<BmRotaryMatrix> rotvec;
  /* read rotary matrix (GCRS to ITRS) from input file */
  rotvec = parse_rotary(argv[2]);

   /* read EOPS (we will need DUT1) */
  //dso::EopSeries eop;
  //{
  //  auto imjd = orbvec[0].epoch.imjd();
  //  const auto t1 = dso::MjdEpoch(imjd - 2);
  //  const auto t2 = dso::MjdEpoch(imjd + 3);
  //  /* parse EOP values to the EopSeries instance for any t >= t1 and t < t2 */
  //  if (dso::parse_iers_C04(argv[3], t1, t2, eop)) {
  //    fprintf(stderr, "ERROR Failed parsing eop file\n");
  //    return 1;
  //  }
  //}

  const auto ieops = parse_eops(argv[3]);

  dso::StokesCoeffs stokes(DEGREE, ORDER);

  /* spit out a title for plotting */
  if (formatD3Plot) {
    printf("mjd,sec,refval,val,component\n");
  } else {
    printf("#title Pole Tide\n");
  }

  /* compare results epoch by epoch */
  //double fargs[6];
  //dso::EopRecord eops;
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 3> g;
  auto acc = accvec.begin();
  auto rot = rotvec.begin();
  auto eop = ieops.begin();
  for (const auto &in : orbvec) {
    /* GPSTime */
    const auto tt = in.epoch;

    /* interpolate EOPs, we need pole coordinates */
    //if (dso::EopSeries::out_of_bounds(eop.interpolate(tt, eops))) {
    //  fprintf(stderr, "Failed to interpolate: Epoch is out of bounds!\n");
    //  return 1;
    //}

    ///* compute fundamental arguments at given epoch */
    //dso::fundarg(tt, fargs);

    //    /* compute gmst using an approximate value for UT1 (linear interpolation) */
    //double dut1_approx;
    //eop.approx_dut1(tt, dut1_approx);
    //const double gmst = dso::gmst(tt, tt.tt2ut1(dut1_approx));

    ///* add libration effect [micro as] */
    //{
    //  double dxp, dyp, dut1, dlod;
    //  dso::deop_libration(fargs, gmst, dxp, dyp, dut1, dlod);
    //  eops.xp() += dxp * 1e-6;
    //  eops.yp() += dyp * 1e-6;
    //}

    ///* add ocean tidal effect [micro as] */
    //{
    //  double dxp, dyp, dut1, dlod;
    //  dso::deop_ocean_tide(fargs, gmst, dxp, dyp, dut1, dlod);
    //  eops.xp() += dxp * 1e-6;
    //  eops.yp() += dyp * 1e-6;
    //}

    /* clear Stokes coefficients */
    stokes.clear();

    /* (xp,yp) rad to arcsec */
    const double xp = dso::rad2sec(eop->xp);
    const double yp = dso::rad2sec(eop->yp);
    assert(eop->epoch == in.epoch);

    /**/
    //dso::poleTide::stokes_coeffs(tt, eops.xp(), eops.yp(), stokes.C(2, 1),
    //                             stokes.S(2, 1));
    dso::poleTide::stokes_coeffs(tt, xp, yp, stokes.C(2, 1), stokes.S(2, 1));

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
