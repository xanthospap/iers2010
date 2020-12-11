#include "iers2010.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>

#define TEST_STATUS_SUCCESS 0
#define TEST_STATUS_FAILURE 1

// default location of file: gpt2_5.grd
const char gpt2grd[] = "/usr/local/share/iers10/gpt2_5.grd";
bool file_exists(const char *str) {
  bool status = true;
  std::ifstream fin(str);
  if (!fin.good())
    status = false;
  fin.close();
  return status;
}

int main() {
  std::cout << "\nDriver program to test the c++ implementation of the IERS2010"
               "\n";

  // testing fundarg
  std::cout << "----------------------------------------\n";
  std::cout << "> fundarg\n";
  std::cout << "----------------------------------------\n";
  double fargs[5];
  iers2010::fundarg(0.230441004960306, fargs);
  for (int i = 0; i < 5; ++i)
    printf("\t%+20.15f\n", fargs[i]);

  // testing pmsdnut2
  std::cout << "----------------------------------------\n";
  std::cout << "> pmsdnut2\n";
  std::cout << "----------------------------------------\n";
  iers2010::pmsdnut2(59961.357706175185740, fargs[0], fargs[1]);
  for (int i = 0; i < 2; ++i)
    printf("\t%+20.15f\n", fargs[i]);

  // testing utlibr
  std::cout << "----------------------------------------\n";
  std::cout << "> utlibr\n";
  std::cout << "----------------------------------------\n";
  iers2010::utlibr(59961.357706175185740, fargs[0], fargs[1]);
  for (int i = 0; i < 2; ++i)
    printf("\t%+20.15f\n", fargs[i]);

  // testing fcnnut
  std::cout << "----------------------------------------\n";
  std::cout << "> fcnnut\n";
  std::cout << "----------------------------------------\n";
  iers2010::fcnnut(59961.357706175185740, fargs[0], fargs[1], fargs[2],
                   fargs[3]);
  for (int i = 0; i < 4; ++i)
    printf("\t%+20.15f\n", fargs[i]);

  // testing arg2
  std::cout << "----------------------------------------\n";
  std::cout << "> arg2\n";
  std::cout << "----------------------------------------\n";
  double arg2_out[11];
  iers2010::arg2(2023, 17.357706175185740, arg2_out);
  for (int i = 0; i < 11; ++i)
    printf("\t%+20.15f\n", arg2_out[i]);

  // testing cnmtx
  std::cout << "----------------------------------------\n";
  std::cout << "> cnmtx\n";
  std::cout << "----------------------------------------\n";
  double cnmtx_out[12];
  iers2010::oeop::cnmtx(59961.357706175185740, cnmtx_out);
  for (int i = 0; i < 12; ++i)
    printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing ortho_eop
  std::cout << "----------------------------------------\n";
  std::cout << "> ortho_eop\n";
  std::cout << "----------------------------------------\n";
  iers2010::ortho_eop(59961.357706175185740, cnmtx_out[0], cnmtx_out[1],
                      cnmtx_out[2]);
  for (int i = 0; i < 3; ++i)
    printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing rg_zont2
  std::cout << "----------------------------------------\n";
  std::cout << "> rg_zont2\n";
  std::cout << "----------------------------------------\n";
  iers2010::rg_zont2(0.230441004960306, cnmtx_out[0], cnmtx_out[1],
                     cnmtx_out[2]);
  for (int i = 0; i < 3; ++i)
    printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing fcul_a
  std::cout << "----------------------------------------\n";
  std::cout << "> fcul_a\n";
  std::cout << "----------------------------------------\n";
  cnmtx_out[0] =
      iers2010::fcul_a(-19.912441285299039, 1210.1513e0, 245.008e0, 26.81e0);
  printf("\t%+20.15f\n", cnmtx_out[0]);

  // testing fcul_b
  std::cout << "----------------------------------------\n";
  std::cout << "> fcul_b\n";
  std::cout << "----------------------------------------\n";
  cnmtx_out[0] = iers2010::fcul_b(-19.912441285299039, 1210.1513e0,
                                  17.357706175185740e0, 26.81e0);
  printf("\t%+20.15f\n", cnmtx_out[0]);

  // testing fculzd_hpa
  std::cout << "----------------------------------------\n";
  std::cout << "> fculzd_hpa\n";
  std::cout << "----------------------------------------\n";
  iers2010::fcul_zd_hpa(-19.912441285299039, 1210.1513e0, 1712.68567e0,
                        20.68274e0, 0.54121e0, cnmtx_out[0], cnmtx_out[1],
                        cnmtx_out[2]);
  for (int i = 0; i < 3; ++i)
    printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing gmf
  std::cout << "----------------------------------------\n";
  std::cout << "> gmf\n";
  std::cout << "----------------------------------------\n";
  iers2010::gmf(59961.357706175185740, -0.347537662538520, 1.648122322251867,
                1210.151e0, 1.1029408012e0, cnmtx_out[0], cnmtx_out[1]);
  for (int i = 0; i < 2; ++i)
    printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing vmf1
  std::cout << "----------------------------------------\n";
  std::cout << "> vmf1\n";
  std::cout << "----------------------------------------\n";
  iers2010::vmf1(0.573496442182711, 0.004487220617989, 59961.357706175185740,
                 -0.347537662538520, 1.1029408012e0, cnmtx_out[0],
                 cnmtx_out[1]);
  for (int i = 0; i < 2; ++i)
    printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing vmf1_ht
  std::cout << "----------------------------------------\n";
  std::cout << "> vmf1_ht\n";
  std::cout << "----------------------------------------\n";
  iers2010::vmf1_ht(0.573496442182711, 0.004487220617989, 59961.357706175185740,
                    -0.347537662538520, 1210.15134e0, 1.1029408012e0,
                    cnmtx_out[0], cnmtx_out[1]);
  for (int i = 0; i < 2; ++i)
    printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing gpt
  std::cout << "----------------------------------------\n";
  std::cout << "> gpt\n";
  std::cout << "----------------------------------------\n";
  iers2010::gpt(59961.357706175185740, -0.347537662538520, 1.648122322251867,
                1210.151e0, cnmtx_out[0], cnmtx_out[1], cnmtx_out[2]);
  for (int i = 0; i < 3; ++i)
    printf("\t%+20.15f\n", cnmtx_out[i]);

  // testing gpt2
  std::cout << "----------------------------------------\n";
  std::cout << "> gpt2\n";
  std::cout << "----------------------------------------\n";
  if (!file_exists(gpt2grd)) {
    std::cerr << "\n[ERROR] Cannot locate/open the grd file \"gpt2_5.grd\"\n";
    return 1;
  } else {
    iers2010::gpt2(59961.357706175185740, -0.347537662538520, 1.648122322251867,
                   1210.151e0, 1, cnmtx_out[0], cnmtx_out[1], cnmtx_out[2],
                   cnmtx_out[3], cnmtx_out[4], cnmtx_out[5], cnmtx_out[6],
                   gpt2grd);
    for (int i = 0; i < 7; ++i)
      printf("\t%+20.15f\n", cnmtx_out[i]);
  }

  return 0;
}
