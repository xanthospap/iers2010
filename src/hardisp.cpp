#include "hardisp.hpp"
#include "ggdatetime/datetime_write.hpp"
#include "iers2010.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>

using iers2010::hisp::nl;
using iers2010::hisp::nt;
using iers2010::hisp::ntin;

/// Argv are:
///                  -# year (integer) (Note 4)
///                  -# [month (integer) day_of_month (integer)], or
///                     day_of_year (integer)
///                  -# hour (integer)
///                  -# minute (integer)
///                  -# second (integer)
///                  -# Number of output epochs to be written out/computed
///                     (integer)
///                  -# Sampling interval (in seconds) (double)
int main(int argc, char *argv[]) {
  //  Check number of arguments from command line, then read them in
  if (argc < 8 || argc > 9) {
    printf(" Usage:\n");
    printf("   hardisp yr [d-of-yr | month day] hr min sec num samp\n");
    printf(" Where \n");
    printf("   the UTC date given is the time of the first term output\n");
    printf("   num is the number of output epochs to be written out\n");
    printf("   samp is the sample interval (seconds)\n");
    printf(" The harmonics file (amp and phase of displacement) is \n");
    printf("   read from standard input in the BLQ format used by  \n");
    printf("   Scherneck and Bos\n");
    printf(" Results are written to standard output (units = m):\n");
    printf("      dU    dS    dW   \n");
    printf("   using format: 3F14.6 \n");
    return 1;
  }

  int next = 1;
  int irnt;
  double samp;
  double it[6];
  try {
    it[0] = std::stoi(argv[next++]);
    if (argc == 7) { /* day of year provided */
      it[1] = std::stoi(argv[next++]);
    } else { /* month - day of month provided */
      it[1] = std::stoi(argv[next++]);
      it[2] = std::stoi(argv[next++]);
    }
    it[3] = std::stoi(argv[next++]);
    it[4] = std::stoi(argv[next++]);
    it[5] = (int)std::stof(argv[next++]);
    irnt = std::stoi(argv[next++]);
    samp = std::stod(argv[next++]);
  } catch (std::invalid_argument &) {
    std::cerr
        << "[ERROR] Invalid argument while reading input arguments. Fatal.\n";
    return 1;
  } catch (std::out_of_range &) {
    std::cerr
        << "[ERROR] Invalid argument while reading input arguments. Fatal.\n";
    return 1;
  }

  /// Construct the date
  dso::hours _h(it[3]);
  dso::minutes _m(it[4]);
  dso::seconds _s(it[5]);
  dso::datetime<dso::seconds> epoch;
  if (argc == 8) {
    epoch = dso::datetime<dso::seconds>{dso::year(it[0]),
                                        dso::day_of_year(it[1]), _h, _m, _s};
  } else {
    epoch = dso::datetime<dso::seconds>{dso::year(it[0]),
                                        dso::month(it[1]),
                                        dso::day_of_month(it[2]),
                                        _h,
                                        _m,
                                        _s};
  }

  /// Read in ocean loading coefficients from stdin
  double tamp[3][ntin], tph[3][ntin];
  if (iers2010::hisp::read_hardisp_args(tamp, tph)) {
    std::cerr << "[ERROR] Failed to read harmonics. Fatal.\n";
    return 1;
  }

  double *du = new double[irnt];
  double *dw = new double[irnt];
  double *ds = new double[irnt];
  iers2010::hisp::hardisp_impl(irnt, samp, tamp, tph, epoch, du, ds, dw);
  for (int i = 0; i < irnt; i++)
    printf("\n%14.6f %14.6f %14.6f", du[i], ds[i], dw[i]);
  printf("\n");
  delete[] du;
  delete[] dw;
  delete[] ds;
}
