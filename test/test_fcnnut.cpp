#include "iers2010.hpp"
#include <cassert>
#include <cmath>

const double PRECISION = 1e-9;

int main() {
  double fcx, fcy, fcdx, fcdy;
    const double imjd[] = {
.547900000000000E+05,
.457000000000000E+05,
.457000000000001E+05,
.457001234500000E+05,
.457009999900000E+05,
.500821234567890E+05,
.500829999999999E+05,
.500830000000009E+05,
.500831234567890E+05,
.537358888888888E+05,
.537359888888888E+05,
.537359000000000E+05,
.537360000000000E+06,
.537360000000001E+05,
.537360000000001E+06,
.537360000000011E+06,
.537360000000111E+06,
.537360000011111E+06,
.541010000000000E+05,
.541020000000000E+05,
.562930000000000E+05,
.562939990000000E+05,
.562939999999999E+05,
.567930000000000E+05   };
iers2010::fcnnut(imjd[  0], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-(-0.176801229006627E+03 ))<PRECISION);
 assert(std::abs(fcy-(-0.935185530890376E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.374557377049180E+01 ))<PRECISION);
assert(std::abs(fcdy-( 0.374557377049180E+01 ))<PRECISION);
iers2010::fcnnut(imjd[  1], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-(-0.225784766311316E+02 ))<PRECISION);
 assert(std::abs(fcy-( 0.291377983557002E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.394400000000000E+02 ))<PRECISION);
assert(std::abs(fcdy-( 0.394400000000000E+02 ))<PRECISION);
iers2010::fcnnut(imjd[  2], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-(-0.225784766310629E+02 ))<PRECISION);
 assert(std::abs(fcy-( 0.291377983557711E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.394399999999952E+02 ))<PRECISION);
assert(std::abs(fcdy-( 0.394399999999952E+02 ))<PRECISION);
iers2010::fcnnut(imjd[  3], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-(-0.224952007104674E+02 ))<PRECISION);
 assert(std::abs(fcy-( 0.292235445826739E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.394341985245902E+02 ))<PRECISION);
assert(std::abs(fcdy-( 0.394341985245902E+02 ))<PRECISION);
iers2010::fcnnut(imjd[  4], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-(-0.218971607531333E+02 ))<PRECISION);
 assert(std::abs(fcy-( 0.298264247846512E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.393930059344264E+02 ))<PRECISION);
assert(std::abs(fcdy-( 0.393930059344264E+02 ))<PRECISION);
iers2010::fcnnut(imjd[  5], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.721303027132816E+02 ))<PRECISION);
 assert(std::abs(fcy-(-0.548478769720419E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.153427857263692E+02 ))<PRECISION);
assert(std::abs(fcdy-( 0.153427857263692E+02 ))<PRECISION);
iers2010::fcnnut(imjd[  6], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.714341558926489E+02 ))<PRECISION);
 assert(std::abs(fcy-(-0.557641809040944E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.153400000000003E+02 ))<PRECISION);
assert(std::abs(fcdy-( 0.153400000000003E+02 ))<PRECISION);
iers2010::fcnnut(imjd[  7], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.714341558917906E+02 ))<PRECISION);
 assert(std::abs(fcy-(-0.557641809052142E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.153399999999839E+02 ))<PRECISION);
assert(std::abs(fcdy-( 0.153399999999839E+02 ))<PRECISION);
iers2010::fcnnut(imjd[  8], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.713277144000154E+02 ))<PRECISION);
 assert(std::abs(fcy-(-0.559028293406938E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.153377939688523E+02 ))<PRECISION);
assert(std::abs(fcdy-( 0.153377939688523E+02 ))<PRECISION);
iers2010::fcnnut(imjd[  9], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.146727979861791E+03 ))<PRECISION);
 assert(std::abs(fcy-(-0.636382477508325E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.438032267884348E+01 ))<PRECISION);
assert(std::abs(fcdy-( 0.438032267884348E+01 ))<PRECISION);
iers2010::fcnnut(imjd[ 10], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.146645529962806E+03 ))<PRECISION);
 assert(std::abs(fcy-(-0.638485946811882E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.438003226788458E+01 ))<PRECISION);
assert(std::abs(fcdy-( 0.438003226788458E+01 ))<PRECISION);
iers2010::fcnnut(imjd[ 11], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.146718833653121E+03 ))<PRECISION);
 assert(std::abs(fcy-(-0.636616248080545E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.438029041095890E+01 ))<PRECISION);
assert(std::abs(fcdy-( 0.438029041095890E+01 ))<PRECISION);
iers2010::fcnnut(imjd[ 12], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.255059481391124E+03 ))<PRECISION);
 assert(std::abs(fcy-(-0.116623356805106E+03 ))<PRECISION);
assert(std::abs(fcdx-( 0.127486475000000E+06 ))<PRECISION);
assert(std::abs(fcdy-( 0.127486475000000E+06 ))<PRECISION);
iers2010::fcnnut(imjd[ 13], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.146636350251330E+03 ))<PRECISION);
 assert(std::abs(fcy-(-0.638719600839845E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.437999999999982E+01 ))<PRECISION);
assert(std::abs(fcdy-( 0.437999999999982E+01 ))<PRECISION);
iers2010::fcnnut(imjd[ 14], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.255059481389320E+03 ))<PRECISION);
 assert(std::abs(fcy-(-0.116623356809049E+03 ))<PRECISION);
assert(std::abs(fcdx-( 0.127486475000000E+06 ))<PRECISION);
assert(std::abs(fcdy-( 0.127486475000000E+06 ))<PRECISION);
iers2010::fcnnut(imjd[ 15], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.255059481372244E+03 ))<PRECISION);
 assert(std::abs(fcy-(-0.116623356846397E+03 ))<PRECISION);
assert(std::abs(fcdx-( 0.127486475000003E+06 ))<PRECISION);
assert(std::abs(fcdy-( 0.127486475000003E+06 ))<PRECISION);
iers2010::fcnnut(imjd[ 16], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.255059481201898E+03 ))<PRECISION);
 assert(std::abs(fcy-(-0.116623357218950E+03 ))<PRECISION);
assert(std::abs(fcdx-( 0.127486475000029E+06 ))<PRECISION);
assert(std::abs(fcdy-( 0.127486475000029E+06 ))<PRECISION);
iers2010::fcnnut(imjd[ 17], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.255059462465773E+03 ))<PRECISION);
 assert(std::abs(fcy-(-0.116623398195520E+03 ))<PRECISION);
assert(std::abs(fcdx-( 0.127486475002944E+06 ))<PRECISION);
assert(std::abs(fcdy-( 0.127486475002944E+06 ))<PRECISION);
iers2010::fcnnut(imjd[ 18], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.129198729667723E+03 ))<PRECISION);
 assert(std::abs(fcy-( 0.111810219355149E+03 ))<PRECISION);
assert(std::abs(fcdx-( 0.374000000000000E+01 ))<PRECISION);
assert(std::abs(fcdy-( 0.374000000000000E+01 ))<PRECISION);
iers2010::fcnnut(imjd[ 19], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.130769015964554E+03 ))<PRECISION);
 assert(std::abs(fcy-( 0.110017195221890E+03 ))<PRECISION);
assert(std::abs(fcdx-( 0.373928767123288E+01 ))<PRECISION);
assert(std::abs(fcdy-( 0.373928767123288E+01 ))<PRECISION);
iers2010::fcnnut(imjd[ 20], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.170350869482061E+03 ))<PRECISION);
 assert(std::abs(fcy-( 0.222793464147192E+03 ))<PRECISION);
assert(std::abs(fcdx-( 0.372000000000000E+01 ))<PRECISION);
assert(std::abs(fcdy-( 0.372000000000000E+01 ))<PRECISION);
iers2010::fcnnut(imjd[ 21], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.173583250765591E+03 ))<PRECISION);
 assert(std::abs(fcy-( 0.220284364977749E+03 ))<PRECISION);
assert(std::abs(fcdx-( 0.398473500000091E+01 ))<PRECISION);
assert(std::abs(fcdy-( 0.398473500000091E+01 ))<PRECISION);
iers2010::fcnnut(imjd[ 22], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.173586467983635E+03 ))<PRECISION);
 assert(std::abs(fcy-( 0.220281829783953E+03 ))<PRECISION);
assert(std::abs(fcdx-( 0.398499999997301E+01 ))<PRECISION);
assert(std::abs(fcdy-( 0.398499999997301E+01 ))<PRECISION);
iers2010::fcnnut(imjd[ 23], fcx, fcy, fcdx, fcdy);
 assert(std::abs(fcx-( 0.279020492410669E+03 ))<PRECISION);
 assert(std::abs(fcy-(-0.283533281098987E+02 ))<PRECISION);
assert(std::abs(fcdx-( 0.136220000000000E+03 ))<PRECISION);
assert(std::abs(fcdy-( 0.136220000000000E+03 ))<PRECISION);
}
