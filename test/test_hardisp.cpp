#include "blqstrm.hpp"
#include "datetime/datetime_write.hpp"
#include "datetime/dtcalendar.hpp"
#include "hardisp.hpp"
#include "iers2010.hpp"
#include <stdio.h>

using iers2010::BlqIn;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "\n[ERROR] Need to provide an input BLQ file";
    return 1;
  }

  // create a BlqIn instance
  BlqIn blq(argv[1]);

  // find and read the station we want in the BLQ file (note that we change
  // phase sign to immidiately use hardisp)
  if (!blq.find_station("ONSA")) {
    std::cerr << "\n[ERROR] Could not find station ONS in the BLQ file";
    return 1;
  }
  double a1[3][11], a2[3][11];
  std::string station;
  if (blq.read_next_station(station, a1, a2, true)) {
    std::cerr << "\n[ERROR] Failed to read records for station \"" << station
              << "\"";
  }

  // call hardisp
  dso::datetime d(dso::year(2009), dso::month(6), dso::day_of_month(25),
                  dso::hours(1), dso::minutes(10), dso::seconds(45));

  // number of samples and sample interval in seconds
  const int irnt = 24;
  const double sample_sec = 3600e0;
  double *displacements = new double[3 * irnt];
  double *dw = displacements, *du = displacements + irnt,
         *ds = displacements + 2 * irnt;
  iers2010::hisp::hardisp_impl(irnt, sample_sec, a1, a2, d, du, ds, dw);
  for (int i = 0; i < irnt; i++)
    printf("\n%14.6f %14.6f %14.6f", du[i], ds[i], dw[i]);
  printf("\n");
  delete[] displacements;

  printf("\n");
  return 0;
}
