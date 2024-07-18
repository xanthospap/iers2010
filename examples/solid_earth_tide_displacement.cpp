#include "fundarg.hpp"
#include "solid_earth_tide.hpp"
#include <charconv>
#include <cstdio>
#include <cstring>
#include "datetime/datetime_read.hpp"

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc != ) {
    fprintf(stderr, "Usage: X Y Z T_START T_END [eopc04.1962-now] [de421.bsp] [naif*.tls]\n");
    return 1;
  }

  /* parse station coordinates */
  Eigen::Matrix<double, 3, 1> rsta;
  {
    double data[3];
    for (int i=0; i<3; i++) {
      auto s = std::strlen(argv[i]);
      auto res = std::from_chars(argv[i], argv[i]+s, data[i]);
      if (res.errc != std::errc{}) {
        fprintf(stderr, "Error! Failed parsing coordinates\n");
        return 2;
      }
    }
    rsta(0) = data[0];
    rsta(1) = data[1];
    rsta(2) = data[2];
  }

  MjdEpoch start, end;
  {
    char buf[128];
    std::strcpy(buf, argv[3]);
    auto sz = std::strlen(buf);
    buf[sz] = ' ';
    buf[sz+1] = '\0';
    std::strcat(buf, argv[4]);
    TwoPartDateUTC ustart = from_utc_char<YMDFormat::YYYYMMDD, HMSFormat::HHMMSSF>(buf);
    
    std::strcpy(buf, argv[5]);
    sz = std::strlen(buf);
    buf[sz] = ' ';
    buf[sz+1] = '\0';
    std::strcat(buf, argv[6]);
    TwoPartDateUTC uend = from_utc_char<YMDFormat::YYYYMMDD, HMSFormat::HHMMSSF>(buf);

    start = ustart.utc2tai();
    end = uend.utc2tai();
  }

  /* EOP info */
  EopSeries eop;
  {
    auto imjd = start.imjd();
    const auto t1 = MjdEpoch(imjd - 2);
    imjd = end.imjd();
    const auto t2 = MjdEpoch(imjd + 3);
    /* parse EOP values to the EopSeries instance for any t >= t1 and t < t2 */
    if (parse_iers_C04(argv[7], t1, t2, eop)) {
      fprintf(stderr, "ERROR Failed parsing eop file\n");
      return 1;
    }
  }

   /* Load CSPICE/NAIF Kernels */
  dso::load_spice_kernel(argv[8]);
  dso::load_spice_kernel(argv[9]);

  dso::CartesianCrd rsta, rmon, rsun;
  SolidEarthTide set(3.986004418e14, 6378136.6e0, 4.9048695e12,
                     1.32712442099e20);
  
  TwoPartDateUTC utc(dso::datetime<nanoseconds>(
      year(2009), month(4), day_of_month(13), nanoseconds(0)));
  MjdEpoch tt(utc.utc2tt());
  MjdEpoch ut(tt.tt2ut1(0.3089055e0));

  /* Fundamental arguments (IERS 2003) */
  double fa[14] = {
      fal03(tt),  falp03(tt), faf03(tt),  fad03(tt),  faom03(tt),
      fame03(tt), fave03(tt), fae03(tt),  fama03(tt), faju03(tt),
      fasa03(tt), faur03(tt), fane03(tt), fapa03(tt),
  };

  auto dr = set.displacement(tt, ut, rsta.mv, rmon.mv, rsun.mv, fa);
