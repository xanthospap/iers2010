#include <stdio.h>
#include "ggdatetime/dtcalendar.hpp"
#include "ggdatetime/datetime_write.hpp"
#include "iers2010.hpp"
#include "blqstrm.hpp"
#include "hardisp.hpp"

using iers2010::BlqIn;

int main (int argc,char *argv[])
{
  if (argc<1) {
    std::cerr<<"\n[ERROR] Need to provide an input BLQ file";
    return 1;
  }

  // create a BlqIn instance
  BlqIn blq (argv[1]);

  // find and read the station we ant in the BLQ file (note that we change
  // phase sign to immidiately use hardisp)
  if (!blq.find_station("ONSA")) {
    std::cerr<<"\n[ERROR] Could not find station ONS in the BLQ file";
    return 1;
  }
  double a1[3][11], a2[3][11];
  std::string station;
  if (blq.read_next_station(station, a1, a2, true)) {
    std::cerr<<"\n[ERROR] Failed to read records for station \""<<station<<"\"";
  }

  // call hardisp
  ngpt::datetime d(ngpt::year(2009), ngpt::month(6), ngpt::day_of_month(25),
    ngpt::hours(1), ngpt::minutes(10), ngpt::seconds(45));
  /*
  std::cout<<"\nDatetime: "<<ngpt::strftime_ymd_hmfs(d);
  std::cout<<"\n24, 3600e0";
  for (int i=0; i<3; i++) {
    std::cout<<"\n";
    for (int j=0; j<11; j++) std::cout<<" "<<a1[i][j];
  }
  for (int i=0; i<3; i++) {
    std::cout<<"\n";
    for (int j=0; j<11; j++) std::cout<<" "<<a2[i][j];
  }
  */
  iers2010::hisp::hardisp_impl(24, 3600e0, a1, a2, d);
  
  printf ("\n");
  return 0;
}
