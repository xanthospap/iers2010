#include "icgemio.hpp"
#include <cstdio>

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [ICGEM file]\n", argv[0]);
    return 1;
  }

  Icgem gfc(argv[1]);
  if (gfc.parse_header()) {
    fprintf(stderr, "Failed to parse ICGEM header for %s\n",
            gfc.filename().c_str());
    return 1;
  }

  printf("Icgem filename: %s\n", gfc.filename().c_str());
  printf("Product Type  : %s\n", gfc.product_type());
  printf("Model Name    : %s\n", gfc.model_name());
  printf("Tide System   : %s\n", gfc.tide_system());
  printf("Errors        : %s\n", gfc.errors());
  printf("GM            : %.3f[m^3s^-2]\n", gfc.gm());
  printf("R_earth       : %.15f[m]\n", gfc.radius());
  printf("Max degree    : %d\n", gfc.max_degree());

  return 0;
}
