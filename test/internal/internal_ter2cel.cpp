#include "cel2ter.hpp"
#include "sofa.h"
#include <charconv>
#include <fstream>
#include <vector>

struct GraceOrbit {
  dso::TwoPartDate t;
  Eigen::Matrix<double, 3, 1> pos;
  Eigen::Matrix<double, 3, 1> vel;
};

namespace {
const char *skip_ws(const char *c) {
  while (*c && *c == ' ')
    ++c;
  return c;
}
} /* namespace */

int estate_line(const char *line, GraceOrbit &orb);
int parse_orbit(const char *fn, std::vector<GraceOrbit> &orb);
Eigen::Matrix<double, 3, 1> sofa_gcrf2itrf(const Eigen::Matrix<double, 3, 1> &r,
                                           const dso::Itrs2Gcrs &R);

int verbose = 0;
int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s [C04/14 EOP file] [ORBIT FILE (ECI)]\n",
            argv[0]);
    return 1;
  }

  /* store orbit here */
  std::vector<GraceOrbit> orb;

  /* parse orbit file */
  if (parse_orbit(argv[2], orb)) {
    fprintf(stderr, "ERROR. Failed parsing orbit file %s\n", argv[2]);
    return 1;
  }

  /* Create an EOP LookUp table */
  dso::EopLookUpTable eop_lut;
  {
    const int start = (int)orb[0].t.big() - 5;
    const int end = (int)orb[0].t.big() + 5;
    if (dso::parse_iers_C04(argv[1], start, end, eop_lut)) {
      fprintf(stderr, "ERROR. Failed collecting EOP data\n");
      return 1;
    }
    eop_lut.utc2tt();
    eop_lut.regularize(false);
  }

  /* an ITRF-to-GCRF Rotation for general use */
  dso::Itrs2Gcrs R(orb[0].t, &eop_lut);

  /* keep results here */
  double max_dx = 0e0;
  double max_dy = 0e0;
  double max_dz = 0e0;
  double max_dr = 0e0;
  double max_dx_sofa = 0e0;
  double max_dy_sofa = 0e0;
  double max_dz_sofa = 0e0;
  double max_dr_sofa = 0e0;

  /* one by one, transform sp3 input orbit from ITRF to GCRF and back
   * check results
   */
  for (const auto &i : orb) {
    R.prepare(i.t);
    /* GCRF-to-ITRF  (i.pos->r1) */
    Eigen::Matrix<double, 3, 1> r1 = R.gcrf2itrf(i.pos);
    /* ITRF-to-GCRF  (r1->r2) */
    Eigen::Matrix<double, 3, 1> r2 = R.itrf2gcrf(r1);

    /* sofa result */
    Eigen::Matrix<double, 3, 1> rs2 = sofa_gcrf2itrf(i.pos, R);

    if (verbose) {
      printf("%.12e %+.9e %+.9e %+.9e %.9e\n", R.tt().as_mjd(),
             i.pos(0) - r2(0), i.pos(1) - r2(1), i.pos(2) - r2(2),
             (i.pos - r2).norm());
    }

    const auto dr = r2 - i.pos;
    if (std::abs(dr(0)) > std::abs(max_dx))
      max_dx = dr(0);
    if (std::abs(dr(1)) > std::abs(max_dy))
      max_dy = dr(1);
    if (std::abs(dr(2)) > std::abs(max_dz))
      max_dz = dr(2);
    if (dr.norm() > max_dr)
      max_dr = dr.norm();
    
    const auto dr_sofa = rs2 - i.pos;
    if (std::abs(dr_sofa(0)) > std::abs(max_dx_sofa))
      max_dx_sofa = dr_sofa(0);
    if (std::abs(dr_sofa(1)) > std::abs(max_dy_sofa))
      max_dy_sofa = dr_sofa(1);
    if (std::abs(dr_sofa(2)) > std::abs(max_dz_sofa))
      max_dz_sofa = dr_sofa(2);
    if (dr_sofa.norm() > max_dr_sofa)
      max_dr_sofa = dr_sofa.norm();
  }

  printf("Function         #Tests #Fails #Maxerror[mm/sec] Status Type\n");
  printf("---------------------------------------------------------------\n");
  printf("%8s %7s %6d %6d %+.1e %s %s\n", "internal_ter2cel", "dX",
         (int)orb.size(), 0, max_dx * 1e3, "-", "linear");
  printf("%8s %7s %6d %6d %+.1e %s %s\n", "internal_ter2cel", "dY",
         (int)orb.size(), 0, max_dy * 1e3, "-", "linear");
  printf("%8s %7s %6d %6d %+.1e %s %s\n", "internal_ter2cel", "dZ",
         (int)orb.size(), 0, max_dz * 1e3, "-", "linear");
  printf("%8s %7s %6d %6d %+.1e %s %s\n", "internal_ter2cel", "ds",
         (int)orb.size(), 0, max_dr * 1e3, "-", "linear");
  printf("%8s %7s %6d %6d %+.1e %s %s\n", "internal_ter2celIau", "dX",
         (int)orb.size(), 0, max_dx_sofa * 1e3, "-", "linear");
  printf("%8s %7s %6d %6d %+.1e %s %s\n", "internal_ter2celIau", "dY",
         (int)orb.size(), 0, max_dy_sofa * 1e3, "-", "linear");
  printf("%8s %7s %6d %6d %+.1e %s %s\n", "internal_ter2celIau", "dZ",
         (int)orb.size(), 0, max_dz_sofa * 1e3, "-", "linear");
  printf("%8s %7s %6d %6d %+.1e %s %s\n", "internal_ter2celIau", "ds",
         (int)orb.size(), 0, max_dr_sofa * 1e3, "-", "linear");

  return 0;
}

int estate_line(const char *line, GraceOrbit &orb) {
  double d[10];
  const char *start = line;
  const char *end = line + std::strlen(line);
  for (int i = 0; i < 10; i++) {
    auto res = std::from_chars(skip_ws(start), end, d[i]);
    start = res.ptr;
    if (res.ec != std::errc{})
      return 1;
  }
  double ip, fp;
  fp = std::modf(d[0], &ip);
  orb.t = dso::TwoPartDate(ip, fp);
  orb.pos(0) = d[1];
  orb.pos(1) = d[2];
  orb.pos(2) = d[3];
  orb.vel(0) = d[4];
  orb.vel(1) = d[5];
  orb.vel(2) = d[6];
  return 0;
}

int parse_orbit(const char *fn, std::vector<GraceOrbit> &orb) {
  orb.clear();

  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR. Failed to find input file %s\n", fn);
    return 1;
  }

  constexpr const int MC = 1024;
  char line[MC];

  while (fin.getline(line, MC)) {
    GraceOrbit o;
    if (!estate_line(line, o))
      orb.emplace_back(o);
  }

  return 0;
}

/* SOFA GCRS to ITRS */
Eigen::Matrix<double, 3, 1> sofa_gcrf2itrf(const Eigen::Matrix<double, 3, 1> &r,
                                           const dso::Itrs2Gcrs &R) {
  /* for SOFA, transform MJD to JD */
  const double tt1 = R.tt().big() + dso::mjd0_jd;
  const double tt2 = R.tt().small();

  /* CIP and CIO, IAU 2006/2000A. */
  double x, y;
  iauXy06(tt1, tt2, &x, &y);
  double s = iauS06(tt1, tt2, x, y);

  /* Add CIP corrections (CIP offsets wrt IAU 2006/2000A in [rad]) */
  x += dso::sec2rad(R.eop().dx);
  y += dso::sec2rad(R.eop().dy);

  /* GCRS to CIRS matrix. */
  double rc2i[3][3];
  iauC2ixys(x, y, s, rc2i);

  /* Earth rotation Angle. */
  const double ut1 = R.ut1().big() + dso::mjd0_jd;
  const double ut2 = R.ut1().small();
  double era = iauEra00(ut1, ut2);

  /* Form celestial-terrestrial matrix (no polar motion yet). */
  double rc2ti[3][3];
  iauCr(rc2i, rc2ti);
  iauRz(era, rc2ti);

  /* Polar motion matrix (TIRS->ITRS, IERS 2003). */
  double rpom[3][3];
  double sp = iauSp00(tt1, tt2);
  iauPom00(R.eop().xp * iers2010::DAS2R, R.eop().yp * iers2010::DAS2R, sp,
           rpom);

  /* Form celestial-terrestrial matrix (including polar motion). */
  double rc2it[3][3];
  iauRxr(rpom, rc2ti, rc2it);

  /* copy to output matrix */
  double Rc2i[3][3];
  iauCr(rc2it, Rc2i);

  /* transform GCRS vector to ITRS vector */
  double ritrf[3];
  double rgcrf[3] = {r(0), r(1), r(2)};
  iauRxp(rc2it, rgcrf, ritrf);

  /* transpose the matrix to get the inverse traansformation */
  double Rt2c[3][3];
  iauTr(Rc2i, Rt2c);
  /* transform ITRS vector to GCRS vector */
  iauRxp(Rt2c, ritrf, rgcrf);

  Eigen::Matrix<double, 3, 1> _r;
  _r << rgcrf[0], rgcrf[1], rgcrf[2];

  return _r;
}
