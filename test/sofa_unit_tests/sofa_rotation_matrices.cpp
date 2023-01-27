#include "iau.hpp"
#include "sofa.hpp"
#include <cstdio>
#include <cassert>
#include "unit_test_help.hpp"

using iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

int main(int argc, char *argv[]) {
    if (argc>1) {
        fprintf(stderr, "Ignoring command line arguments!\n");
    }

    double am, as;
    double rot[3][3];
    Eigen::Matrix<double, 3, 3> r1,r2,r3,r4,r5;
    double m1[3][3], m2[3][3], m3[3][3], m4[3][3], m5[3][3];
    double dm[20],ds[20];
    for (int i=0; i<NUM_TESTS; i++) {
        // random date (MJD, TT)
        const auto tt = random_mjd();
        const auto jc = tt.jcenturies_sinceJ2000();

        // add a few seconds to TT to get a random UT1 date
        const auto ut = add_random_seconds(tt, -60e0, 60e0);

        // split TT JD in a "J2000" fashion
        auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();
        
        iauPnm00a(jdtt._big,jdtt._small,rot);
        assert(approx_equal(pnm00a(tt), rot));

        const double dpsi = random_angle(-dso::DPI/10, dso::DPI/10);
        const double deps = random_angle(-dso::DPI/10, dso::DPI/10);
        pn00a(tt,dpsi,deps,am,r1,r2,r3,r4,r5);
        iauPn00a(jdtt._big,jdtt._small,dpsi,deps,as,r1,r2,r3,r4,r5);
        assert(approx_equal(r1,m1));
        assert(approx_equal(r2,m2));
        assert(approx_equal(r3,m3));
        assert(approx_equal(r4,m4));
        assert(approx_equal(r5,m5));
        assert(approx_equal(am,as));
        
        pn00(tt,dpsi,deps,am,r1,r2,r3,r4,r5);
        iauPn00(jdtt._big,jdtt._small,dpsi,deps,as,r1,r2,r3,r4,r5);
        assert(approx_equal(r1,m1));
        assert(approx_equal(r2,m2));
        assert(approx_equal(r3,m3));
        assert(approx_equal(r4,m4));
        assert(approx_equal(r5,m5));
        assert(approx_equal(am,as));
        
        pn06(tt,dpsi,deps,am,r1,r2,r3,r4,r5);
        iauPn06(jdtt._big,jdtt._small,dpsi,deps,as,r1,r2,r3,r4,r5);
        assert(approx_equal(r1,m1));
        assert(approx_equal(r2,m2));
        assert(approx_equal(r3,m3));
        assert(approx_equal(r4,m4));
        assert(approx_equal(r5,m5));
        assert(approx_equal(am,as));
        
        p06e(tt,dm[0],dm[1],dm[2],dm[3],dm[4],dm[5],dm[6],dm[7],dm[8],dm[9],dm[10],dm[11],dm[12],dm[13],dm[14],dm[15]);
        iauP06e(jdtt._big,jdtt._small,ds[0],ds[1],ds[2],ds[3],ds[4],ds[5],ds[6],ds[7],ds[8],ds[9],ds[10],ds[11],ds[12],ds[13],ds[14],ds[15]);
        for (int i=0;i<16;i++) assert(approx_equal(dm[i],ds[i]));
        
        bp00(tt,r1,r2,r3);
        iauBp00(jdtt._big,jdtt._small,m1,m2,m3);
        assert(approx_equal(r1,m1));
        assert(approx_equal(r2,m2));
        assert(approx_equal(r3,m3));
    }
}