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

/* results and input from ARG2;
 * last update: 2011 October 07
 */
double arg2_res[] = {
    2.849663065753787805e0,
    6.28318080000000023e0,
    4.926040134021299366e0,
    1.608450491115348768e0,
    2.375021572352622456e0,
    0.4746414933980958040e0,
    3.908159227647345801e0,
    2.551018561669245344e0,
    5.041990012540757959e0,
    4.206816878908014701e0,
    1.608463638294885811e0
};

/* results from CNMTX;
 * last update: 2010 March 17
 */
double cnmtx_res[] = {
    15.35873641938967360e0,
    9.784941251812741214e0,
    -5.520740128266865554e0,
    3.575314211234633888e0,
    -13.93717453496387648e0,
    -9.167400321705855504e0,
    5.532815475865292321e0,
    9.558741883500834646e0,
    -10.22541212627272600e0,
    0.8367570529461261231e0,
    1.946355176475630611e0,
    -13.55702062247304696e0
};

/* results from ORTHO_EOP;
 * last update 2010 March 19
 */
double ortho_eop_res[] = {
    -162.8386373279636530e0,
    117.7907525842668974e0,
    -23.39092370609808214e0
};

/* results from RG_ZONT2;
 * last update 2011 December 20
 */
double rg_zont2_res[] = {
    7.983287678576557467E-002,
    5.035331113978199288E-005,
    -4.249711616463017E-014
};

int main()
{
    int status          = TEST_STATUS_SUCCESS;
    int fundarg_status  = TEST_STATUS_SUCCESS;
    int pmsdnut2_status = TEST_STATUS_SUCCESS;
    int utlibr_status   = TEST_STATUS_SUCCESS;
    int fcnnut_status   = TEST_STATUS_SUCCESS;
    int arg2_status     = TEST_STATUS_SUCCESS;
    int cnmtx_status    = TEST_STATUS_SUCCESS;
    int ortho_eop_status= TEST_STATUS_SUCCESS;
    int rg_zont2_status = TEST_STATUS_SUCCESS;
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
    
    // testing arg2
    std::cout<<"----------------------------------------\n";
    std::cout<<"> arg2\n";
    std::cout<<"----------------------------------------\n";
    double arg2_array[11];
    iers2010::arg2(2008, 311.5e0, arg2_array);
    for (int i=0; i<11; ++i) {
        diff = std::fabs(arg2_array[i]-arg2_res[i]); 
        std::cout << "\tdiff: " << diff <<"\n";
        if (diff > 1e-15) {
            arg2_status = TEST_STATUS_FAILURE;
        }
    }
    std::cout << "status: " << (arg2_status ? "failed!\n" : "ok\n" );
    status += arg2_status;

    // testing cnmtx
    std::cout<<"----------------------------------------\n";
    std::cout<<"> cnmtx\n";
    std::cout<<"----------------------------------------\n";
    double cnmtx_tmp[12];
    iers2010::oeop::cnmtx(54964.0e0, cnmtx_tmp);
    for (int i=0; i<12; ++i) {
        diff = std::fabs(cnmtx_tmp[i]-cnmtx_res[i]);
        std::cout << "\targument[" << i << "] diff: " << diff <<"\n";
        if (diff > 1e-15) {
            cnmtx_status = TEST_STATUS_FAILURE;
        }
    }
    std::cout << "status: " << (cnmtx_status ? "failed!\n" : "ok\n" );
    status += cnmtx_status;
    
    // testing ortho_eop
    std::cout<<"----------------------------------------\n";
    std::cout<<"> ortho_eop\n";
    std::cout<<"----------------------------------------\n";
    iers2010::ortho_eop(47100.0e0, cnmtx_tmp[0], cnmtx_tmp[1], cnmtx_tmp[2]);
    for (int i=0; i<3; ++i) {
        diff = std::fabs(cnmtx_tmp[i]-ortho_eop_res[i]);
        std::cout << "\targument[" << i << "] diff: " << diff <<"\n";
        if (diff > 1e-15) {
            ortho_eop_status = TEST_STATUS_FAILURE;
        }
    }
    std::cout << "status: " << (ortho_eop_status ? "failed!\n" : "ok\n" );
    status += ortho_eop_status;
    
    // testing rg_zont2
    std::cout<<"----------------------------------------\n";
    std::cout<<"> rg_zont2\n";
    std::cout<<"----------------------------------------\n";
    iers2010::rg_zont2(.07995893223819302e0, cnmtx_tmp[0], cnmtx_tmp[1], cnmtx_tmp[2]);
    for (int i=0; i<3; ++i) {
        diff = std::fabs(cnmtx_tmp[i]-rg_zont2_res[i]);
        std::cout << "\targument[" << i << "] diff: " << diff <<"\n";
        if (diff > 1e-15) {
            rg_zont2_status = TEST_STATUS_FAILURE;
        }
    }
    std::cout << "status: " << (rg_zont2_status ? "failed!\n" : "ok\n" );
    status += rg_zont2_status;
    
    return status;    
}
