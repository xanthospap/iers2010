#include "earth_rotation.hpp"
#include "geodesy/units.hpp"
#include "iau.hpp"
#include <array>
#include <cstdio>
#include <random>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <fstream>
#include <charconv>
#include <cstring>

constexpr const double PTOLERANCE = 9e-2;
constexpr const double VTOLERANCE = 9e-5;

const char *mskipws(const char *line) noexcept {
  while (line && *line == ' ')
    ++line;
  return line;
}

struct RefOrbit {
  dso::MjdEpoch t;
  /* terrestrial */
  double xi, yi, zi;
  double vxi, vyi, vzi;
  /* celestial */
  double xc, yc, zc;
  double vxc, vyc, vzc;

  RefOrbit(const dso::MjdEpoch &t_, const double* const d_) :
  t(t_), xi(d_[0]), yi(d_[1]), zi(d_[2]), vxi(d_[3]), vyi(d_[4]), vzi(d_[5]),
  xc(d_[6]), yc(d_[7]), zc(d_[8]), vxc(d_[9]), vyc(d_[10]), vzc(d_[11]){}
}; /* RefOrbit */

std::vector<RefOrbit> parse_ref_orbit(const char *fn) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR. Failed openning input file %s\n", fn);
    return std::vector<RefOrbit>{};
  }

  std::vector<RefOrbit> vec;
  double data[13];
  char line[512];
  while (fin.getline(line, 512)) {
    int error = 0;
    int sz = std::strlen(line);
    const char *str = line;
    for (int i = 0; i < 13; i++) {
      auto res = std::from_chars(mskipws(str), line + sz, data[i]);
      if (res.ec != std::errc{})
        ++error;
      str = res.ptr;
    }
    if (!error) {
      int imjd = (int)data[0];
      double fsec = (data[0] - imjd) * 86400e0;
      const dso::MjdEpoch t(imjd, dso::FractionalSeconds{fsec});
      vec.emplace_back(t, data+1);
    }
  }
  
  return vec;
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s [ref transformation data file] [EOP]\n", argv[0]);
    return 1;
  }

  /* parse data */
  const auto ref = parse_ref_orbit(argv[1]);
  if (!ref.size()) return 9;

  /* first and last epochs for eops */
  const auto t0 = ref[0].t.add_seconds(dso::FractionalSeconds(-2*86400));
  const auto t1 = (ref.end()-1)->t.add_seconds(dso::FractionalSeconds(2*86400));

  /* create an instance to hold EOP series */
  dso::EopSeries eops;

  /* parse EOP values to the EopSeries instance for any t >= t1 and t < t2 */
  if (dso::parse_iers_C04(argv[2], t0, t1, eops)) {
    fprintf(stderr, "ERROR Failed parsing eop file\n");
    return 1;
  }

  /* regularize the series (UT1 and LOD) */
  eops.regularize();

  double fargs[14];
  for (const auto &d : ref) {
    /* interpolate EOPS */
    dso::EopRecord eopr;
    if (dso::EopSeries::out_of_bounds(
            eops.interpolate(d.t.tai2tt(), eopr))) {
      fprintf(stderr, "Failed to interpolate: Epoch is out of bounds!\n");
      return 1;
    }
    
    /* compute (X,Y) CIP and fundamental arguments */
    double Xcip, Ycip;
    dso::xycip06a(d.t.tai2tt(), Xcip, Ycip, fargs);

    /* use fundamental arguments to compute s */
    const double s = dso::s06(d.t.tai2tt(), Xcip, Ycip, fargs);

    /* apply CIP corrections */
    Xcip += dso::sec2rad(eopr.dX());
    Ycip += dso::sec2rad(eopr.dY());

    /* compute gmst using an approximate value for UT1 (linear interpolation) */
    double dut1_approx;
    eops.approx_dut1(d.t.tai2tt(), dut1_approx);
    const double gmst =
        dso::gmst(d.t.tai2tt(), d.t.tai2ut1(dut1_approx));

    /* add libration effect [micro as] */
    {
      double dxp, dyp, dut1, dlod;
      dso::deop_libration(fargs, gmst, dxp, dyp, dut1, dlod);
      eopr.xp() += dxp * 1e-6;
      eopr.yp() += dyp * 1e-6;
      eopr.dut() += dut1 * 1e-6;
      eopr.lod() += dlod * 1e-6;
    }

    /* add ocean tidal effect [micro as] */
    {
      double dxp, dyp, dut1, dlod;
      dso::deop_ocean_tide(fargs, gmst, dxp, dyp, dut1, dlod);
      eopr.xp() += dxp * 1e-6;
      eopr.yp() += dyp * 1e-6;
      eopr.dut() += dut1 * 1e-6;
      eopr.lod() += dlod * 1e-6;
    }

    /* de-regularize */
    {
      double ut1_cor;
      double lod_cor;
      double omega_cor;

      dso::deop_zonal_tide(fargs, ut1_cor, lod_cor, omega_cor);
      /* apply (note: microseconds to seconds) */
      eopr.dut() += (ut1_cor * 1e-6);
      eopr.lod() += (lod_cor * 1e-6);
    }

    Eigen::Vector3d ri, vi;
    ri << d.xi, d.yi, d.zi;
    vi << d.vxi, d.vyi, d.vzi;
    /* transform position (ITRS-to-GCRS)  - Option 1 */
    const auto q1 = dso::c2i06a(d.t.tai2tt(), eopr);
    Eigen::Vector3d rc = q1.conjugate() * ri;
    // printf("%+.3f %+.3f %+.3f [P1]\n", rc(0)-d.xc, rc(1)-d.yc, rc(2)-d.zc);
    assert(std::abs(rc(0)-d.xc)<PTOLERANCE);
    assert(std::abs(rc(1)-d.yc)<PTOLERANCE);
    assert(std::abs(rc(2)-d.zc)<PTOLERANCE);

    /* transform position (ITRS-to-GCRS)  - Option 2 */
    const auto q2 = dso::c2i06a_bz(d.t.tai2tt(), eopr);
    rc = q2.conjugate() * ri;
    // printf("%+.3f %+.3f %+.3f [P2]\n", rc(0)-d.xc, rc(1)-d.yc, rc(2)-d.zc);
    assert(std::abs(rc(0)-d.xc)<PTOLERANCE);
    assert(std::abs(rc(1)-d.yc)<PTOLERANCE);
    assert(std::abs(rc(2)-d.zc)<PTOLERANCE);

    /* transform position (ITRS-to-GCRS)  - Option 3 */
    Eigen::Quaterniond q_c2tirs, q_tirs2i;
    const auto omega = dso::c2i06a(d.t.tai2tt(), eopr, q_c2tirs, q_tirs2i);
    rc = q_c2tirs.conjugate() * (q_tirs2i.conjugate() * ri);
    assert(std::abs(rc(0)-d.xc)<PTOLERANCE);
    assert(std::abs(rc(1)-d.yc)<PTOLERANCE);
    assert(std::abs(rc(2)-d.zc)<PTOLERANCE);
    /* transform velocity (ITRS-to-GCRS)  - Option 3 */
    auto vc = q_c2tirs.conjugate() * (q_tirs2i.conjugate() * vi + omega.cross(q_tirs2i.conjugate() * ri));
    assert(std::abs(vc(0)-d.vxc)<VTOLERANCE);
    assert(std::abs(vc(1)-d.vyc)<VTOLERANCE);
    assert(std::abs(vc(2)-d.vzc)<VTOLERANCE);
    //printf("%+.3f %+.3f %+.3f [P3]\n", rc(0)-d.xc, rc(1)-d.yc, rc(2)-d.zc);
    //printf("%+.3f %+.3f %+.3f [V3]\n", vc(0)-d.vxc, vc(1)-d.vyc, vc(2)-d.vzc);
    
    /* transform position (ITRS-to-GCRS) - Option 4 */
    Eigen::Matrix<double,3,3> dMdt;
    const auto q3 = dso::c2i06a(d.t.tai2tt(), eopr, dMdt);
    rc = q3.conjugate() * ri;
    assert(std::abs(rc(0)-d.xc)<PTOLERANCE);
    assert(std::abs(rc(1)-d.yc)<PTOLERANCE);
    assert(std::abs(rc(2)-d.zc)<PTOLERANCE);
    /* transform velocity (ITRS-to-GCRS) - Option 4 */
    vc = q3.conjugate() * vi + dMdt.transpose() * ri;
    assert(std::abs(vc(0)-d.vxc)<VTOLERANCE);
    assert(std::abs(vc(1)-d.vyc)<VTOLERANCE);
    assert(std::abs(vc(2)-d.vzc)<VTOLERANCE);
    //printf("%+.3f %+.3f %+.3f [P4]\n", rc(0)-d.xc, rc(1)-d.yc, rc(2)-d.zc);
    //printf("%+.3f %+.3f %+.3f [V4]\n", vc(0)-d.vxc, vc(1)-d.vyc, vc(2)-d.vzc);
  }

  return 0;
}

