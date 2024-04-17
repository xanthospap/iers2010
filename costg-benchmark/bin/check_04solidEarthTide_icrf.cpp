#include "costg_utils.hpp"
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "gravity.hpp"
#include "icgemio.hpp"
#include "planets.hpp"
#include "iersconst.hpp"
#include "solid_earth_tide.hpp"
#include "eop.hpp"
#include "fundarg.hpp"

using namespace costg;
constexpr const double GM_Sun = 1.32712442076e20;
constexpr const double GM_Moon = 0.49028010560e13;
constexpr const int DEGREE = 4;

int main(int argc, char *argv[]) {
  if (argc != 7) {
    fprintf(
        stderr,
        "Usage: %s [00orbit_icrf.txt] [04solidEarthTide_icrf.txt] [de421.bsp] [naif*.tls] [eopc04.1962-now] [01earthRotation_rotaryMatrix.txt]\n",
        argv[0]);
    return 1;
  }

  /* Load CSPICE/NAIF Kernels */
  dso::load_spice_kernel(argv[3]);
  dso::load_spice_kernel(argv[4]);

  /* read orbit from input file */
  const auto orbvec = parse_orbit(argv[1]);
  printf("#Note: read %d data sets from input file\n", (int)orbvec.size());

  /* read acceleration from input file */
  const auto accvec = parse_acceleration(argv[2]);
  printf("#Note: read %d data sets from input file\n", (int)accvec.size());
  
  /* read rotary matrix (GCRS to ITRS) from input file */
  const auto rotvec = parse_rotary(argv[6]);
  printf("#Note: read %d data sets from input file\n", (int)accvec.size());

  /* read and parse EOP iformation */
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
  dso::SolidEarthTide setide(iers2010::GMe, iers2010::Re, GM_Moon, GM_Sun);

    /* allocate scratch space for computations */
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> W(DEGREE + 3,
                                                                    DEGREE + 3);
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> M(DEGREE + 3,
                                                                    DEGREE + 3);

  /* compare results epoch by epoch */
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 1> rsun,rmon,rsat;
  Eigen::Matrix<double, 3, 3> g;
  double fargs[5];
  auto acc = accvec.begin();
  auto rot = rotvec.begin();
  for (const auto &in : orbvec) {
    /* GPSTime to TT */
    const auto t = in.epoch.gps2tai().tai2tt();

    /* rotation matrix */
    const Eigen::Matrix<double,3,3> R = rot->R;
    assert(rot->epoch == in.epoch);
    
    /* get Sun+Moon position in ICRF */
    if (dso::planet_pos(dso::Planet::SUN, t, rsun)) {
      fprintf(stderr, "ERROR Failed to compute Sun position!\n");
      return 2;
    }
    if (dso::planet_pos(dso::Planet::MOON, t, rmon)) {
      fprintf(stderr, "ERROR Failed to compute Sun position!\n");
      return 2;
    }

    /* GCRF to ITRF */
    rsun = R * rsun;
    rmon = R * rmon;
    rsat = R * in.xyz;
    
    /* compute gmst using an approximate value for UT1 (linear interpolation) */
    double dut1_approx;
    eop.approx_dut1(t, dut1_approx);

    /* compute fundamental (Delaunay) arguments fot t */
    dso::fundarg(t, fargs);
    
    /* compute solid earth tide geopotential variations */
    setide.stokes_coeffs(t, t.tt2ut1(dut1_approx), rmon, rsun, fargs);

    /* compute acceleration for given epoch/position */
    dso::sh2gradient_cunningham(setide.stokes_coeffs(), rsat, a, g, -1, -1, -1,
                                -1, &W, &M);

    /* acceleration: ITRF to GCRF */
    a = R.transpose() * a;

    /* get COSTG result */
    if (acc->epoch != in.epoch) {
      fprintf(stderr, "ERROR Failed to match epochs in input files\n");
      return 1;
    }

    printf("%d %.9f %.15e %.15e %.15e %.15e %.15e %.15e\n", in.epoch.imjd(),
           in.epoch.seconds(), acc->axyz(0), acc->axyz(1), acc->axyz(2), a(0),
           a(1), a(2));
    ++acc;
    ++rot;
  }
}
