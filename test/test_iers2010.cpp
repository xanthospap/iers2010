#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include "iers2010.hpp"

#define TEST_STATUS_SUCCESS 0
#define TEST_STATUS_FAILURE 1


/* results for GPT.F;
 * last update 2011 October 18
 */
double gpt_res[] =
{
  918.0710638757363995e0,
  19.31914181012882992e0,
  -42.19185643717770517e0
};

/* results for GPT2
 * last update 31 May 2013
 */
double gpt2_res_a[] =
{
  1002.555031899637015,
  22.121274188593230,
  -6.525262725051438,
  15.625300656392879,
  0.001264668783781,
  0.000572557214095,
  44.05609719704228
};

double gpt2_res_b[] =
{
  1003.487051814098550,
  11.952642050310883,
  -5.469287999554646,
  9.578015648003742,
  0.001239543515998,
  0.000556030804004,
  44.056097197042284
};

// default location of file: gpt2_5.grd
const char gpt2grd[] = "/usr/local/share/iers10/gpt2_5.grd";
bool
file_exists(const char* str)
{
  bool status = true;
  std::ifstream fin(str);
  if (!fin.good()) status = false;
  fin.close();
  return status;
}

int main()
{
  int status          = TEST_STATUS_SUCCESS;
  int gpt_status      = TEST_STATUS_SUCCESS;
  int gpt2_status     = TEST_STATUS_SUCCESS;
  double diff;

  // format results
  std::cout.setf(std::ios::fixed/*, std:: ios::floatfield*/);
  std::cout.precision(15);
  std::cerr.setf(std::ios::fixed/*, std:: ios::floatfield*/);
  std::cerr.precision(15);

  std::cout<<"\nDriver program to test the c++ implementation of the IERS2010"
    "\n";

  // testing fundarg
  std::cout<<"----------------------------------------\n";
  std::cout<<"> fundarg\n";
  std::cout<<"----------------------------------------\n";
  double fargs[5];
  iers2010::fundarg(   0.097212536121941, fargs);
  for (int i=0; i<5; ++i) printf("\t%+20.15f\n", fargs[i]);

  // testing pmsdnut2
  std::cout<<"----------------------------------------\n";
  std::cout<<"> pmsdnut2\n";
  std::cout<<"----------------------------------------\n";
  iers2010::pmsdnut2(55095.187881853897125, fargs[0], fargs[1]);
  for (int i=0; i<2; ++i) printf("\t%+20.15f\n", fargs[i]);

  // testing utlibr
  std::cout<<"----------------------------------------\n";
  std::cout<<"> utlibr\n";
  std::cout<<"----------------------------------------\n";
  iers2010::utlibr(55095.187881853897125, fargs[0], fargs[1]);
  for (int i=0; i<2; ++i) printf("\t%+20.15f\n", fargs[i]);

  // testing fcnnut
  std::cout<<"----------------------------------------\n";
  std::cout<<"> fcnnut\n";
  std::cout<<"----------------------------------------\n";
  iers2010::fcnnut(55095.187881853897125, fargs[0], fargs[1], fargs[2], fargs[3]);
  for (int i=0; i<4; ++i) printf("\t%+20.15f\n", fargs[i]);

  // testing arg2
  std::cout<<"----------------------------------------\n";
  std::cout<<"> arg2\n";
  std::cout<<"----------------------------------------\n";
  double arg2_out[11];
  iers2010::arg2(2009,  264.187881853897125, arg2_out);
  for (int i=0; i<11; ++i) printf("\t%+20.15f\n", arg2_out[i]);

  // testing cnmtx
  std::cout<<"----------------------------------------\n";
  std::cout<<"> cnmtx\n";
  std::cout<<"----------------------------------------\n";
  double cnmtx_out[12];
  iers2010::oeop::cnmtx(55095.187881853897125, cnmtx_out);
  for (int i=0; i<12; ++i) printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing ortho_eop
  std::cout<<"----------------------------------------\n";
  std::cout<<"> ortho_eop\n";
  std::cout<<"----------------------------------------\n";
  iers2010::ortho_eop(55095.187881853897125, cnmtx_out[0], cnmtx_out[1], cnmtx_out[2]);
  for (int i=0; i<3; ++i) printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing rg_zont2
  std::cout<<"----------------------------------------\n";
  std::cout<<"> rg_zont2\n";
  std::cout<<"----------------------------------------\n";
  iers2010::rg_zont2(   0.097212536121941, cnmtx_out[0], cnmtx_out[1], cnmtx_out[2]);
  for (int i=0; i<3; ++i) printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing fcul_a
  std::cout<<"----------------------------------------\n";
  std::cout<<"> fcul_a\n";
  std::cout<<"----------------------------------------\n";
  cnmtx_out[0] = iers2010::fcul_a(  30.576398284960717,  -636.0803e0,    300.918e0, 78.04e0);
  printf("\t%+20.15f\n", cnmtx_out[0]);

  // testing fcul_b
  std::cout<<"----------------------------------------\n";
  std::cout<<"> fcul_b\n";
  std::cout<<"----------------------------------------\n";
  cnmtx_out[0] = iers2010::fcul_b(  30.576398284960717,  -636.0803e0,  264.187881853897125e0, 78.04e0);
  printf("\t%+20.15f\n", cnmtx_out[0]);

  // testing fculzd_hpa
  std::cout<<"----------------------------------------\n";
  std::cout<<"> fculzd_hpa\n";
  std::cout<<"----------------------------------------\n";
  iers2010::fcul_zd_hpa(  30.576398284960717,  -636.0803e0,  636.90271e0,   28.46009e0,    0.54769e0,
      cnmtx_out[0], cnmtx_out[1], cnmtx_out[2]);
  for (int i=0; i<3; ++i) printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing gmf
  std::cout<<"----------------------------------------\n";
  std::cout<<"> gmf\n";
  std::cout<<"----------------------------------------\n";
  iers2010::gmf(55095.187881853897125,    0.533658823473712,   -0.416959733920350,   -636.080e0,    0.2087427636e0,
      cnmtx_out[0], cnmtx_out[1]);
  for (int i=0; i<2; ++i) printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing vmf1
  std::cout<<"----------------------------------------\n";
  std::cout<<"> vmf1\n";
  std::cout<<"----------------------------------------\n";
  iers2010::vmf1(   0.145915391735760,    0.413168631696726, 55095.187881853897125,    0.533658823473712,    0.2087427636e0,
      cnmtx_out[0], cnmtx_out[1]);
  for (int i=0; i<2; ++i) printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing vmf1_ht
  std::cout<<"----------------------------------------\n";
  std::cout<<"> vmf1_ht\n";
  std::cout<<"----------------------------------------\n";
  iers2010::vmf1_ht(   0.145915391735760,    0.413168631696726, 55095.187881853897125,    0.533658823473712,      -636.08031e0,    0.2087427636e0,
      cnmtx_out[0], cnmtx_out[1]);
  for (int i=0; i<2; ++i) printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing gpt
  std::cout<<"----------------------------------------\n";
  std::cout<<"> gpt\n";
  std::cout<<"----------------------------------------\n";
  iers2010::gpt(55055e0, 0.6708665767e0, -1.393397187e0, 812.546e0,
      cnmtx_out[0], cnmtx_out[1], cnmtx_out[2]);
  for (int i=0; i<3; ++i) {
    diff = std::fabs(cnmtx_out[i]-gpt_res[i]);
    std::cout << "\targument[" << i << "] diff: " << diff <<"\n";
    if (diff > 1e-15) {
      gpt_status = TEST_STATUS_FAILURE;
    }
  }
  std::cout << "status: " << (gpt_status ? "failed!\n" : "ok\n" );
  status += gpt_status;

  // testing gpt2
  double gpt_input[] = { 0.8412486994612668e0,0.28571039855147173e0, 156.e0 };
  std::cout<<"----------------------------------------\n";
  std::cout<<"> gpt2 (Test Case A)\n";
  std::cout<<"----------------------------------------\n";
  if (!file_exists(gpt2grd)) {
    std::cerr<<"\n[ERROR] Cannot locate/open the grd file \"gpt2_5.grd\"\n";
    gpt2_status=1;
  } else {
    int gpt2stat = iers2010::gpt2(56141.e0, gpt_input, gpt_input+1,gpt_input+2, 1, 0,
        cnmtx_out, cnmtx_out+1, cnmtx_out+2, cnmtx_out+3, cnmtx_out+4,
        cnmtx_out+5, cnmtx_out+6, gpt2grd);
    if  (gpt2stat != 0) {
      std::cerr<<"\ngpt2 (case A) failed with status=" << gpt2stat;
      gpt2_status = TEST_STATUS_FAILURE;
    } else {
      for (int i=0; i<7; ++i) {
        diff = std::fabs(cnmtx_out[i]-gpt2_res_a[i]);
        std::cout << "\targument[" << i << "] diff: " << diff <<"\n";
        if (diff > 1e-15) {
          gpt2_status = TEST_STATUS_FAILURE;
        }
      }
    }
  }
  std::cout << "status: " << (gpt2_status ? "failed!\n" : "ok\n" );
  status += gpt2_status;

  // testing gpt2
  std::cout<<"----------------------------------------\n";
  std::cout<<"> gpt2 (Test Case B)\n";
  std::cout<<"----------------------------------------\n";
  if (!file_exists(gpt2grd)) {
    std::cerr<<"\n[ERROR] Cannot locate/open the grd file \"gpt2_5.grd\"\n";
    gpt2_status=1;
  } else {
    int gpt2stat = iers2010::gpt2(56141.e0, gpt_input, gpt_input+1,gpt_input+2, 1, 1,
        cnmtx_out, cnmtx_out+1, cnmtx_out+2, cnmtx_out+3, cnmtx_out+4,
        cnmtx_out+5, cnmtx_out+6, gpt2grd);
    if (gpt2stat != 0) {
      std::cerr<<"\ngpt2 (case B) failed with status=" << gpt2stat;
      gpt2_status = TEST_STATUS_FAILURE;
    } else {
      for (int i=0; i<7; ++i) {
        diff = std::fabs(cnmtx_out[i]-gpt2_res_b[i]);
        std::cout << "\targument[" << i << "] diff: " << diff <<"\n";
        if (diff > 1e-15) {
          gpt2_status = TEST_STATUS_FAILURE;
        }
      }
    }
  }
  std::cout << "status: " << (gpt2_status ? "failed!\n" : "ok\n" );
  status += gpt2_status;

  return status;    
}
