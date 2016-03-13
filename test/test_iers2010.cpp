#include <iostream>
#include <iomanip>
#include "iers2010.hpp"

#define TEST_STATUS_SUCCESS 0
#define TEST_STATUS_FAILURE 1

/* results and input from FUNDARG.F;
 * last update: 2010 February 25
 */
double fundarg_inp { 
    0.07995893223819302e0 };
double fundarg_res[] = {
    2.291187512612069099e0,
    6.212931111003726414e0,
    3.658025792050572989e0,
    4.554139562402433228e0,
    -0.5167379217231804489e0
};

/* results and input from PMSDNUT2;
 * last update: 2011 October  13
 */
double pmsdnut2_inp {
    54335e0 };
double pmsdnut2_res[] = {
    24.83144238273364834e0,
    -14.09240692041837661e0
};

/* results and input from UTLIBR;
 * last update: 2010 June  23
 * (two test-cases)
 */
double utlibr_inp[] = {
    44239.1,           /* test case A */
    55227.4 };         /* test case B */
double utlibr_res[] = {
    2.441143834386761746e0,  /* test case A */
    -14.78971247349449492e0,
    -2.655705844335680244e0, /* test case B */
    27.39445826599846967e0
};

/* results and input from FCNNUT;
 * last update: 2013 December 19
 */
double fcnnut_inp = {
    54790e0 };
double fcnnut_res[] = {
    -176.8012290066270680e0,
    -93.51855308903756736e0,
    3.745573770491803067e0,
    3.745573770491803067e0
};

int main()
{
    int status          = TEST_STATUS_SUCCESS;
    int fundarg_status  = TEST_STATUS_SUCCESS;
    int pmsdnut2_status = TEST_STATUS_SUCCESS;
    int utlibr_status   = TEST_STATUS_SUCCESS;
    int fcnnut_status   = TEST_STATUS_SUCCESS;
    double diff;
    
    // format results
    std::cout.setf(std::ios::fixed/*, std:: ios::floatfield*/);
    std::cout.precision(15);
    std::cerr.setf(std::ios::fixed/*, std:: ios::floatfield*/);
    std::cerr.precision(15);

    // testing fundarg
    std::cout<<"----------------------------------------\n";
    std::cout<<"> fundarg\n";
    std::cout<<"----------------------------------------\n";
    double fargs[5];
    iers2010::fundarg(fundarg_inp, fargs);
    for (int i=0; i<5; ++i) {
        diff = std::fabs(fargs[i]-fundarg_res[i]);
        std::cout << "\targument[" << i << "] diff: " << diff <<"\n";
        if (diff > 1e-15) {
            fundarg_status = TEST_STATUS_FAILURE;
        }
    }
    std::cout << "status: " << (fundarg_status ? "failed!\n" : "ok\n" );
    status += fundarg_status;
    
    // testing pmsdnut2
    std::cout<<"----------------------------------------\n";
    std::cout<<"> pmsdnut2\n";
    std::cout<<"----------------------------------------\n";
    iers2010::pmsdnut2(pmsdnut2_inp, fargs[0], fargs[1]);
    for (int i=0; i<2; ++i) {
        diff = std::fabs(fargs[i]-pmsdnut2_res[i]);
        std::cout << "\tdiff: " << diff <<"\n";
        if ( diff > 1e-15) {
            pmsdnut2_status = TEST_STATUS_FAILURE;
        }
    }
    std::cout << "status: " << (pmsdnut2_status ? "failed!\n" : "ok\n" );
    status += pmsdnut2_status;
    
    // testing utlibr
    std::cout<<"----------------------------------------\n";
    std::cout<<"> utlibr\n";
    std::cout<<"----------------------------------------\n";
    // test-case A
    iers2010::utlibr(utlibr_inp[0], fargs[0], fargs[1]);
    for (int i=0; i<2; ++i) {
        diff = std::fabs(fargs[i]-utlibr_res[i]);
        std::cout << "\tdiff: " << diff <<"\n";
        if (diff > 1e-15) {
            utlibr_status = TEST_STATUS_FAILURE;
        }
    }
    // test-case B
    iers2010::utlibr(utlibr_inp[1], fargs[0], fargs[1]);
    for (int i=0; i<2; ++i) {
        diff = std::fabs(fargs[i]-utlibr_res[i+2]); 
        std::cout << "\tdiff: " << diff <<"\n";
        if ( diff > 1e-15) {
            utlibr_status = TEST_STATUS_FAILURE;
        }
    }
    std::cout << "status: " << (utlibr_status ? "failed!\n" : "ok\n" );
    status += utlibr_status;
    
    // testing fcnnut
    std::cout<<"----------------------------------------\n";
    std::cout<<"> fcnnut\n";
    std::cout<<"----------------------------------------\n";
    iers2010::fcnnut(fcnnut_inp, fargs[0], fargs[1], fargs[2], fargs[3]);
    for (int i=0; i<4; ++i) {
        diff = std::fabs(fargs[i]-fcnnut_res[i]); 
        std::cout << "\tdiff: " << diff <<"\n";
        if (diff > 1e-15) {
            fcnnut_status = TEST_STATUS_FAILURE;
        }
    }
    std::cout << "status: " << (fcnnut_status ? "failed!\n" : "ok\n" );
    status += fcnnut_status;
    
    return status;    
}
