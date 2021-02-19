#include "blqstrm.hpp"
#include "iers2010.hpp"
#include <iostream>

/// A simple program to test the implementation of BlqIn class
/// that is the class tha manipulates BLQ format files.

using iers2010::BlqIn;

int main(int argc, char *argv[]) {

  if (argc != 2) {
    std::cerr << "\n[ERROR] Need to provide a valid .BLQ file (argc #" << argc
              << ")";
    return 1;
  }
  BlqIn blq(argv[1]);

  int status = 0;
  std::string staname;
  double a1[3][11], a2[3][11];

  // Read all stations; if we find station "VALA" then read it else ignore
  // amp and phase values
  while (!status) {

    // Peak next station in file
    if ((status = blq.peak_next_station(staname))) {
      if (status > 0)
        std::cerr << "\n[DEBUG] peak_next_station error";
      break;
    }
    std::cout << "\nNext station in BLQ file is: " << staname;

    // if it is VALA, read and save its values (phase/amp)
    if (staname == "VALA") {
      std::string vala;
      status = blq.read_next_station(vala, a1, a2);
      if (status) {
        std::cerr << "\n[DEBUG] read_next_station error";
      }
      std::cout << "\nRead station VALA (" << vala << ")";
      std::cout << "\nHere are the values:";
      for (int i = 0; i < 3; i++) {
        std::cout << "\n";
        for (int j = 0; j < 11; j++)
          std::cout << " " << a1[i][j];
      }
      for (int i = 0; i < 3; i++) {
        std::cout << "\n";
        for (int j = 0; j < 11; j++)
          std::cout << " " << a2[i][j];
      }
    }

    // else read and ignore the values
    if ((status = blq.skip_next_station())) {
      if (status > 0)
        std::cerr << "\n[DEBUG] skip_next_station failed!";
      break;
    }
  }

  // Now here is a simpler way; suppose we want the values for station "DYNG"
  // Lets see if we can find the station in the file
  if (blq.find_station("DYNG")) {
    // we found it so lets collect its records
    std::string dyng;
    status = blq.read_next_station(dyng, a1, a2);
    if (status) {
      std::cerr << "\n[DEBUG] read_next_station error";
    }
    std::cout << "\nRead station DYNG (" << dyng << ")";
    std::cout << "\nHere are the values:";
    for (int i = 0; i < 3; i++) {
      std::cout << "\n";
      for (int j = 0; j < 11; j++)
        std::cout << " " << a1[i][j];
    }
    for (int i = 0; i < 3; i++) {
      std::cout << "\n";
      for (int j = 0; j < 11; j++)
        std::cout << " " << a2[i][j];
    }
  } else {
    std::cerr << "\n[DEBUG] Cannot find station DYNG";
  }

  // Lets see if a non-existent station is found in the file, e.g. "XXXX"
  if (blq.find_station("XXXX")) {
    std::cerr << "\n[DEBUG] It seems we found station \"XXXX\" in the file!";
  } else {
    std::cout << "\n[DEBUG] OK! station \"XXXX\" was not found in the file";
  }

  std::cout << "\n";
  return 0;
}
