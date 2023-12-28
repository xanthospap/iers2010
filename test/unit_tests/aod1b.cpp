#include "aod1b.hpp"
#include "datetime/datetime_write.hpp"
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

  char buf[64];
  /* iterator for the ATM coefficients */
  dso::Aod1bDataBlockIterator<dso::AOD1BCoefficientType::ATM> it(aod);
  it.set_begin();
  int j = 0;
  while (!j) {
    printf("New data block at %s\n",
          dso::to_char<dso::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSSF,
                         dso::nanoseconds>(it.header().mepoch, buf));
    it.skip();
    j = it.advance();
  }

  return 0;
}
