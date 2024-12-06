#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "eop.hpp"
#include "fundarg.hpp"
#include "gravity.hpp"
#include "icgemio.hpp"
#include "iersconst.hpp"
#include "planets.hpp"
#include "solid_earth_tide.hpp"

using namespace costg;
constexpr const double GM_Sun = 1.32712442076e20;
constexpr const double GM_Moon = 0.49028010560e13;
constexpr const int DEGREE = 4;
constexpr const int formatD3Plot = 0;

int main(int argc, char *argv[]) {
  if (argc != 7) {
    fprintf(
        stderr,
        "Usage: %s [00orbit_itrf.txt] [04solidEarthTide_icrf.txt] [de421.bsp] "
        "[naif*.tls] [eopc04.1962-now] [01earthRotation_rotaryMatrix.txt]\n",
        argv[0]);
    return 1;
  }

  /* Load CSPICE/NAIF Kernels */
  dso::load_spice_kernel(argv[3]);
  dso::load_spice_kernel(argv[4]);

  /* read orbit from input file */
  const auto orbvec = parse_orbit(argv[1]);

  /* read acceleration from input file */
  const auto accvec = parse_acceleration(argv[2]);

  /* read rotary matrix (GCRS to ITRS) from input file */
  const auto rotvec = parse_rotary(argv[6]);

  /* read and parse EOP iformation; we will need UT1-UTC info latter on, to
   * compute e.g. GMST
   */
  dso::EopSeries eop;
  {
    auto imjd = orbvec[0].epoch.imjd();
    const auto t1 = dso::MjdEpoch(imjd - 2);
    const auto t2 = dso::MjdEpoch(imjd + 3);
    /* parse EOP values to the EopSeries instance for any t >= t1 and t < t2 */
    if (dso::parse_iers_C04(argv[5], t1, t2, eop)) {
      fprintf(stderr, "ERROR Failed parsing eop file\n");
      return 1;
    }
  }

  /* A solidEarthTide instance */
  dso::SolidEarthTide setide(iers2010::GMe, iers2010::Re, GM_Sun, GM_Moon);

  /* allocate scratch space for computations (optionally, we could let
   * sh2gradient_cunningham allocate the memmory it needs, but allocating only
   * once is better)
   */
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> W(DEGREE + 3,
                                                                    DEGREE + 3);
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> M(DEGREE + 3,
                                                                    DEGREE + 3);

  /* spit out a title for plotting */
  if (formatD3Plot) {
    printf("mjd,sec,refval,val,component\n");
  } else {
    printf("#title Solid Earth Tide (ephermeris: %s)\n", basename(argv[3]));
  }

  /* compare results epoch by epoch */
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 1> rsun, rmon, rsat;
  Eigen::Matrix<double, 3, 3> g;
  double fargs[5];
  auto acc = accvec.begin();
  auto rot = rotvec.begin();
  for (const auto &in : orbvec) {
    /* we will consider GPSTime */
    const auto t = in.epoch;

    /* rotation matrix */
    const Eigen::Matrix<double, 3, 3> R = rot->R;
    assert(rot->epoch == in.epoch);

    /* get Sun+Moon position in ICRF, using (GPST to) TT */
    if (dso::planet_pos(dso::Planet::SUN, t.gps2tai().tai2tt(), rsun)) {
      fprintf(stderr, "ERROR Failed to compute Sun position!\n");
      return 2;
    }
    if (dso::planet_pos(dso::Planet::MOON, t.gps2tai().tai2tt(), rmon)) {
      fprintf(stderr, "ERROR Failed to compute Sun position!\n");
      return 2;
    }

    /* GCRF to ITRF */
    rsun = R * rsun;
    rmon = R * rmon;
    rsat = in.xyz;

    /* compute gmst using an approximate value for UT1 (linear interpolation) */
    double dut1_approx;
    eop.approx_dut1(t.gps2tai().tai2tt(), dut1_approx);

    /* compute fundamental (Delaunay) arguments for t (using TT) */
    dso::fundarg(t.gps2tai().tai2tt(), fargs);

    /* compute solid earth tide geopotential variations (again, using GPST
     * instead of TT).
     */
    setide.stokes_coeffs(t.gps2tai().tai2tt(), t.gps2tai().tai2ut1(dut1_approx),
                         rmon, rsun, fargs);
    // const auto tutc_ = t.gps2utc();
    // dso::MjdEpoch tutc(tutc_.imjd(), tutc_.seconds());
    // setide.stokes_coeffs(t, tutc, rmon, rsun, fargs);

    /* compute acceleration for given epoch/position (ITRF) */
    dso::sh2gradient_cunningham(setide.stokes_coeffs(), rsat, a, g, -1, -1, -1,
                                -1, &W, &M);

    /* acceleration: ITRF to GCRF */
    a = R.transpose() * a;

    /* get COSTG result */
    if (acc->epoch != in.epoch) {
      fprintf(stderr, "ERROR Failed to match epochs in input files\n");
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
}
