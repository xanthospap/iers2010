#include "sofa.h"
#include "geodesy/units.hpp"
#include "iau.hpp"
#include "eigen3/Eigen/Eigen"
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <random>
#include <cstdio>

void compare_pom(double xp, double yp, double sp) {
  double rpom[3][3];
  iauPom00(xp,yp,sp,rpom);
}

Eigen::Matrix<double, 3, 3> fromq(double era, double s, double sp, double Xcip, double Ycip, double xp, double yp) {
  const double tw = std::cos(xp/2e0) * std::cos(yp/2e0);
  const double aw = std::cos(xp/2e0) * std::sin(yp/2e0);
  const double bw = std::cos(yp/2e0) * std::sin(xp/2e0);
  const double cw = std::sin(xp/2e0) * std::sin(yp/2e0);

  const double thetap = era + sp - s;
  const double ct2 = std::cos(thetap/2);
  const double st2 = std::sin(thetap/2);

/* Obtain the spherical angles E and d. */
  const double r2 = Xcip*Xcip + Ycip*Ycip;
  const double e = (r2 > 0.0) ? std::atan2(Ycip, Xcip) : 0.0;
  const double d = std::atan(std::sqrt(r2 / (1.0 - r2)));
  const double X = std::sin(d) * std::cos(e);
  const double Y = std::sin(d) * std::sin(e);
  const double Z = std::cos(d);

  const double fac = 1e0/(std::sqrt(2e0*(1e0+Z)));
  const double T = fac*(ct2*(X*bw - Y*aw + (1e0+Z)*tw) + st2*(X*aw+Y*bw+(1e0+Z)*cw));
  const double A = fac*(ct2*(X*cw + Y*tw + (1e0+Z)*aw) + st2*(-X*tw+Y*cw-(1e0+Z)*bw));
  const double B = fac*(ct2*(-X*tw + Y*cw + (1e0+Z)*bw) + st2*(-X*cw-Y*tw+(1e0+Z)*aw));
  const double C = fac*(ct2*(-X*aw - Y*bw + (1e0+Z)*cw) + st2*(X*bw-Y*aw-(1e0+Z)*tw));

  Eigen::Quaterniond q(T, A, B, C);
  return q.toRotationMatrix();
}

/* SOFA GCRS to ITRS */
Eigen::Matrix<double, 3, 3> sofa(double jd1, double jd2, double xp, double yp) {

  /* CIP and CIO, IAU 2006/2000A. */
  double x, y;
  iauXy06(jd1, jd2, &x, &y);
  double s = iauS06(jd1, jd2, x, y);

  /* GCRS to CIRS matrix. */
  double rc2i[3][3];
  iauC2ixys(x, y, s, rc2i);

  /* Earth rotation Angle. */
  double era = iauEra00(jd1, jd2);

  /* Form celestial-terrestrial matrix (no polar motion yet). */
  double rc2ti[3][3];
  iauCr(rc2i, rc2ti);
  iauRz(era, rc2ti);

  /* Polar motion matrix (TIRS->ITRS, IERS 2003). */
  double rpom[3][3];
  double sp = iauSp00(jd1, jd2);
  iauPom00(xp, yp, sp, rpom);

  /* Form celestial-terrestrial matrix (including polar motion). */
  double rc2it[3][3];
  iauRxr(rpom, rc2ti, rc2it);

  printf("------------------------------------------------------quaternion\n");
  const auto R2 = fromq(era, s, iauSp00(jd1,jd2), x, y, xp, yp);
    for (int i=0; i<3; i++) {
      for (int j=0;j<3;j++) {
        printf(" %+18.15f ", R2(i,j));
      }
      printf("\n");
    }

/* copy to output matrix */
  Eigen::Matrix<double, 3, 3> R;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      R(i,j) = rc2it[i][j];
    }
  }

  return R;
}

int main() {
    std::uniform_real_distribution<double> unif(-dso::DPI/6, dso::DPI/6);
    std::default_random_engine re;
    
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random(40587, 66154);
    const double jd1 = mjd.imjd() + dso::MJD0_JD;
    const double jd2 = mjd.fractional_days();
    const double xp = unif(re);
    const double yp = unif(re);

    const auto Rsofa = sofa(jd1,jd2,xp,yp);
  printf("------------------------------------------------------euler-----\n");
    for (int i=0; i<3; i++) {
      for (int j=0;j<3;j++) {
        printf(" %+18.15f ", Rsofa(i,j));
      }
      printf("\n");
    }

    return 0;
}
