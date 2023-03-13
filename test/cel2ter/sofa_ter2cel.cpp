#include "cel2ter.hpp"
#include "sofa.h"
#include <cassert>

/*
 * TODO transform this to a test
 * Currently, it basically prints a lot of stuf, for examining shit
 */

struct Sp3Orb {
  dso::datetime<dso::nanoseconds> t_;
  Eigen::Matrix<double,3,1> pos_,vel_;
  dso::TwoPartDate t() const noexcept {return dso::TwoPartDate(t_);}
  Eigen::Matrix<double,3,1> pos() const noexcept {return pos_;}
  Eigen::Matrix<double,3,1> vel() const noexcept {return vel_;}
};

const Sp3Orb Data[] = {
    {dso::datetime<dso::nanoseconds>(dso::year(2020), dso::month(12),
                                     dso::day_of_month(27), dso::hours(8),
                                     dso::minutes(1), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({3487.293316, -896.788983, -6824.852203}),
     Eigen::Matrix<double, 3, 1>({-15363.302305, 65807.746967, -16487.845173})}
};

int main(int argc, char *argv[]) {
  if (argc !=2 ) {
    fprintf(stderr, "Usage: %s [C04/14 EOP file] \n", argv[0]);
    return 1;
  }

  const auto t0 = Data[0].t();
  printf("t0 = %.6f MJD\n", t0.mjd());

  // Create an EOP LookUp table
  dso::EopLookUpTable eop_lut;
  {
    const int start = (int)t0.mjd() - 4;
    const int end = (int)t0.mjd() + 5;
    // parse C04 EOPs
    if (dso::parse_iers_C0414(argv[1], start, end, eop_lut)) {
      fprintf(stderr, "ERROR. Failed collecting EOP data\n");
      return 1;
    }
    eop_lut.utc2tt();
    eop_lut.regularize(false);
  }

  // one by one, transform sp3 input orbit from ITRF to GCRF and back
  // check results
  dso::Itrs2Gcrs Rot(t0, &eop_lut);

  // the following should work
  {
    const auto RR = Rot.gcrf2tirs()*Rot.rpom() ;
    printf("--------------------------------------------------------------\n");
    printf("Rot.rpom() * Rot.gcrf2tirs() - Rot.gcrf2itrf()\n");
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        printf("%+.12f ", RR(i, j) - Rot.gcrf2itrf()(i, j));
      }
      printf("\n");
    }
  }
  
  // construct the ITRF-to-GCRF matrix using SOFA
  double rc2i[3][3];
  double rt2c[3][3];
  {
    auto jd = t0.jd_split<dso::TwoPartDate::JdSplitMethod::DT>();
    const auto Eop = Rot.eop();
    double X,Y;
    iauXy06(jd._big, jd._small,&X,&Y);
    double s = iauS06(jd._big, jd._small,X,Y);
    X += dso::sec2rad(Eop.dx);
    Y += dso::sec2rad(Eop.dy);
    iauC2ixys(X,Y,s,rc2i);
    auto jdut1 = Rot.ut1().jd_split<dso::TwoPartDate::JdSplitMethod::DT>();
    double era = iauEra00(jdut1._big, jdut1._small);
    iauRz(era,rc2i);
    double rpom[3][3];
    s = iauSp00(jd._big, jd._small);
    iauPom00(dso::sec2rad(Eop.xp), dso::sec2rad(Eop.yp), s, rpom);
    iauRxr(rc2i,rpom,rt2c);
    printf("--------------------------------------------------------------\n");
    printf("Here is the Rpom matrix according to SOFA\n");
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        printf("%+14.12f ", rpom[i][j]);
      }
      printf("\n");
    }
    printf("Here is the Rpom matrix according to mine\n");
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        printf("%+14.12f ", Rot.rpom()(i, j));
      }
      printf("\n");
    }
    
    printf("--------------------------------------------------------------\n");
    printf("Here is the Rc2i matrix according to SOFA\n");
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        printf("%+14.12f ", rc2i[i][j]);
      }
      printf("\n");
    }
    printf("Here is the Rc2i matrix according to mine\n");
    const Eigen::Matrix<double, 3, 3> MMR =
        Eigen::AngleAxisd(Rot.earth_rotation_angle(),
                          -Eigen::Vector3d::UnitZ()) *
        Rot.rc2i();
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        printf("%+14.12f ", MMR(i, j));
      }
      printf("\n");
    }

    printf("--------------------------------------------------------------\n");
    printf("Here is the matrix according to SOFA\n");
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        printf("%+14.12f ", rt2c[i][j]);
      }
      printf("\n");
    }
    const auto Rmine = Rot.gcrf2itrf();
    printf("Here is the matrix according to mine\n");
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        printf("%+14.12f ", Rmine(i, j));
      }
      printf("\n");
    }

    printf("--------------------------------------------------------------\n");
    printf("Here is the inverse matrix according to SOFA\n");
    iauTr(rt2c, rc2i);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        printf("%+14.12f ", rc2i[i][j]);
      }
      printf("\n");
    }
    const auto RmineInv = Rot.gcrf2itrf().transpose();
    printf("Here is the inverse matrix according to mine\n");
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        printf("%+14.12f ", RmineInv(i, j));
      }
      printf("\n");
    }
  }
  double rp[3],rc[3],rp2[3];
  rp[0] = Data[0].pos()(0);
  rp[1] = Data[0].pos()(1);
  rp[2] = Data[0].pos()(2);
  iauRxp(rc2i, rp, rc);

  const Eigen::Matrix<double, 3, 1> r_gcrf = Rot.itrf2gcrf(Data[0].pos());
  printf("--------------------------------------------------------------\n");
  printf("GCRF [MINE]: %+.12f %+.12f %+.12f\n", r_gcrf(0), r_gcrf(1), r_gcrf(2));
  printf("GCRF [SOFA]: %+.12f %+.12f %+.12f\n", rc[0], rc[1], rc[2]);

  iauRxp(rt2c, rc, rp2);
  const Eigen::Matrix<double, 3, 1> r_itrf = Rot.gcrf2itrf(r_gcrf);
  printf("--------------------------------------------------------------\n");
  printf("ITRF [MINE]: %+.12f %+.12f %+.12f\n", r_itrf(0), r_itrf(1), r_itrf(2));
  printf("ITRF [SOFA]: %+.12f %+.12f %+.12f\n", rp2[0], rp2[1], rp2[2]);

  printf("--------------------------------------------------------------\n");
  printf("ITRF Diffs: %+.9f %+.9f %+.9f\n", Data[0].pos()(0) - r_itrf(0),
         Data[0].pos()(1) - r_itrf(1), Data[0].pos()(2) - r_itrf(2));
  printf("ITRF Diffs: %+.9f %+.9f %+.9f\n", rp[0] - rp2[0], rp[1] - rp2[1],
         rp[2] - rp2[2]);


  // lets also transform the velocity vectors
  Eigen::Matrix<double, 6, 1> y_in;
  y_in.block<3,1>(0,0) = Data[0].pos();
  y_in.block<3,1>(3,0) = Data[0].vel();
  const Eigen::Matrix<double, 6, 1> y_gcrf = Rot.itrf2gcrf(y_in);
  const Eigen::Matrix<double, 6, 1> y_itrf = Rot.gcrf2itrf(y_gcrf);
  printf("--------------------------------------------------------------\n");
  printf("ITRF [SP3 ]: %+.12f %+.12f %+.12f %+.12f %+.12f %+.12f\n", y_in(0),
         y_in(1), y_in(2), y_in(3), y_in(4), y_in(5));
  printf("ITRF [MINE]: %+.12f %+.12f %+.12f %+.12f %+.12f %+.12f\n", y_itrf(0),
         y_itrf(1), y_itrf(2), y_itrf(3), y_itrf(4), y_itrf(5));

  return 0;
}
