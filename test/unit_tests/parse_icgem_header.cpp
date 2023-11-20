#include "icgemio.hpp"
#include <cstdio>

using namespace dso;

char *tostr(Icgem::ErrorModel e, char *buf) {
  switch (e) {
  case Icgem::ErrorModel::Calibrated:
    return std::strcpy(buf, "calibrated");
  case Icgem::ErrorModel::Formal:
    return std::strcpy(buf, "formal");
  case Icgem::ErrorModel::CalibratedAndFormal:
    return std::strcpy(buf, "calibrated_and_formal");
  case Icgem::ErrorModel::No:
    return std::strcpy(buf, "no");
  }
  return nullptr;
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [ICGEM file]\n", argv[0]);
    return 1;
  }

  char buf[24];

  Icgem gfc(argv[1]);
  printf("Icgem filename: %s\n", gfc.filename().c_str());
  printf("Product Type  : %s\n", gfc.product_type());
  printf("Model Name    : %s\n", gfc.model_name());
  printf("Tide System   : %s\n", gfc.tide_system());
  printf("Errors        : %s\n", tostr(gfc.errors(), buf));
  printf("GM            : %.3f[m^3s^-2]\n", gfc.gm());
  printf("R_earth       : %.15f[m]\n", gfc.radius());
  printf("Max degree    : %d\n", gfc.max_degree());

  return 0;
}
