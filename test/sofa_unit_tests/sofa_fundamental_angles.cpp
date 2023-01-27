#include "iau.hpp"
#include "sofa.hpp"
#include <cstdio>
#include <cassert>
#include "unit_test_help.hpp"

using iers2010::sofa;

constexpr const int NUM_TESTS = 100000;

const char* funcs[] = {
    "fapa03",
    "fane03",
    "faur03",
    "fasa03",
    "faju03",
    "fama03",
    "fae03",
    "fave03",
    "fame03",
    "faom03",
    "fad03",
    "faf03",
    "falp03",
    "fal03",
    "ee06a",
    "gmst06",
    "gmst00"
};

int main(int argc, char *argv[]) {
    if (argc>1) {
        fprintf(stderr, "Ignoring command line arguments!\n");
    }

    double am, as;
    for (int i=0; i<NUM_TESTS; i++) {
        // random date (MJD, TT)
        const auto tt = random_mjd();
        const auto jc = tt.jcenturies_sinceJ2000();

        // add a few seconds to TT to get a random UT1 date
        const auto ut = add_random_seconds(tt, -60e0, 60e0);

        am = fapa03(tt);
        as = iauFapa03(jc);
        assert(approx_equal(am,as));

        am = fane03(tt);
        as = iauFane03(jc);
        assert(approx_equal(am,as));

        am = faur03(tt);
        as = iauFaur03(jc);
        assert(approx_equal(am,as));

        am = fasa03(tt);
        as = iauFasa03(jc);
        assert(approx_equal(am,as));

        am = faju03(tt);
        as = iauFaju03(jc);
        assert(approx_equal(am,as));

        am = fama03(tt);
        as = iauFama03(jc);
        assert(approx_equal(am,as));

        am = fae03(tt);
        as = iauFae03(jc);
        assert(approx_equal(am,as));

        am = fave03(tt);
        as = iauFave03(jc);
        assert(approx_equal(am,as));

        am = fame03(tt);
        as = iauFame03(jc);
        assert(approx_equal(am,as));

        am = faom03(tt);
        as = iauFaom03(jc);
        assert(approx_equal(am,as));

        am = fad03(tt);
        as = iauFad03(jc);
        assert(approx_equal(am,as));

        am = faf03(tt);
        as = iauFaf03(jc);
        assert(approx_equal(am,as));

        am = falp03(tt);
        as = iauFalp03(jc);
        assert(approx_equal(am,as));

        am = fal03(tt);
        as = iauFal03(jc);
        assert(approx_equal(am,as));

        am = ee06a(tt);
        as = iauEe06a(jc);
        assert(approx_equal(am,as));

        am = gmst06(ut,tt);
        auto jdtt = tt.jd_split<dso::TwoPartDate::JdSplitMethod::J2000>();
        auto jdut = ut.jd_split<dso::TwoPartDate::JdSplitMethod::DT>();
        as = iauGmst06(jdut._big, jdut._small, jdtt._big, jdtt._small);
        assert(approx_equal(am,as));
        
        am = gmst00(ut,tt);
        as = iauGmst00(jdut._big, jdut._small, jdtt._big, jdtt._small);
        assert(approx_equal(am,as));
    }

    printf("All functions ok!\n");
}