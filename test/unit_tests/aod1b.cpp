#include "aod1b.hpp"
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [AOD1B]\n", argv[0]);
    return 1;
  }

  dso::Aod1bIn aod(argv[1]);
  if (aod.read_header()) {
    fprintf(stderr, "ERROR Failed parsing AOD1B header!\n");
    return 1;
  }

  printf("Agency    : %s\n", aod.agency());
  printf("File Type : %d\n", aod.file_type());
  printf("File Format: %d\n", aod.file_format());
  printf("Number of header records: %d\n", aod.num_header_records());
  printf("Satellite name: %s\n", aod.satellite());
  printf("Sensor: %s\n", aod.sensor());
  printf("Number of data records %d\n", aod.num_data_records());
  printf("Process Level: %s\n", aod.process_level());
  printf("Pressure type: %s\n", aod.pressure_type());
  printf("Max degree : %d\n", aod.max_degree());
  printf("Coeff errors: %d\n", aod.coeff_errors());
  printf("Coeff normalized: %d\n", aod.coeff_normalized());
  printf("C GM: %.15e\n", aod.GM());

  return 0;
}
