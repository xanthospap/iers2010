#include <iostream>
#include "iers2010/iers2010.hpp"
#include "iers2010/blqstrm.hpp"

/// A simple program to test the implementation of BlqIn class
/// that is the class tha manipulates BLQ format files.

using iers2010::BlqIn;

int main()
{
  BlqIn blq("/home/xanthos/Software/iers2010/data/NTUA.BLQ");
  
  int status=0;
  std::string staname;
  double a1[3][11], a2[3][11];

  while (!status) {
    
    if ((status=blq.peak_next_station(staname))) {
      if (status>0) std::cerr<<"\n[DEBUG] peak_next_station error";
      break;
    }
    std::cout<<"\nNext station in BLQ file is: "<<staname;
    
    if (staname=="VALA") {
      std::string vala;
      status = blq.read_next_station(vala, a1, a2);
      if (status) {
        std::cerr<<"\n[DEBUG] read_next_station error";
      }
      std::cout<<"\nRead station VALA ("<<vala<<")";
      std::cout<<"\nHere are the values:";
      for (int i=0; i<3; i++) {
        std::cout<<"\n";
        for (int j=0 ;j<11; j++) std::cout<<" "<<a1[i][j];
      }
      for (int i=0; i<3; i++) {
        std::cout<<"\n";
        for (int j=0 ;j<11; j++) std::cout<<" "<<a2[i][j];
      }
    }

    if ((status=blq.skip_next_station())) {
      if (status>0) std::cerr<<"\n[DEBUG] skip_next_station failed!";
      break;
    }
  }

  return 0;
}
