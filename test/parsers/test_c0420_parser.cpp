#include "eop.hpp"
#include <cstdio>

struct {
  double mjd, xp, yp, dut1, dx, dy, xrt, yrt, lod;
} InEop[] = {
    // MJD         x(")        y(")       UT1-UTC(s)   dX(")      dY(") xrt(")
    // yrt(")     LOD(s)
    {59917.00, 0.135772, 0.189689, -0.0208103, 0.000324, -0.000149, -0.002950,
     0.000075, 0.0002986},
    {59918.00, 0.132596, 0.189020, -0.0210147, 0.000353, -0.000149, -0.002923,
     -0.000693, 0.0001178},
    {59919.00, 0.130378, 0.188733, -0.0210481, 0.000374, -0.000140, -0.002176,
     -0.000055, -0.0000693},
    {59920.00, 0.127809, 0.188951, -0.0208680, 0.000349, -0.000117, -0.002517,
     0.000100, -0.0002855},
    {59921.00, 0.125572, 0.188560, -0.0204864, 0.000307, -0.000091, -0.002262,
     -0.000222, -0.0004652}};

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [EOP 20/C04 File]\n", argv[0]);
    return 1;
  }
  
  dso::EopLookUpTable eop_lut;

  // parse C04 EOPs
  if (dso::parse_iers_C0420(argv[1], 59917, 59922, eop_lut)) {
    fprintf(stderr, "ERROR. Failed collecting EOP data\n");
    return 1;
  }

  assert(eop_lut.size() == 5);
  for (int i=0; i<5; i++) {
    const auto eop = eop_lut.at(i);
    /*
    printf("%.9f - %.9f = %.12e\n", eop.mjd.mjd(), InEop[i].mjd, eop.mjd.mjd() - InEop[i].mjd);
    printf("%.9f - %.9f = %.12e\n", eop.xp,  InEop[i].xp,  eop.xp  - InEop[i].xp       );
    printf("%.9f - %.9f = %.12e\n", eop.yp,  InEop[i].yp,  eop.yp  - InEop[i].yp       );
    printf("%.9f - %.9f = %.12e\n", eop.dut, InEop[i].dut1,eop.dut - InEop[i].dut1     );
    printf("%.9f - %.9f = %.12e\n", eop.lod, InEop[i].lod, eop.lod - InEop[i].lod      );
    printf("%.9f - %.9f = %.12e\n", eop.dx,  InEop[i].dx,  eop.dx  - InEop[i].dx       );
    printf("%.9f - %.9f = %.12e\n", eop.dy,  InEop[i].dy,  eop.dy  - InEop[i].dy       );
    printf("%.9f - %.9f = %.12e\n", eop.xrt, InEop[i].xrt, eop.xrt - InEop[i].xrt      );
    printf("%.9f - %.9f = %.12e\n", eop.yrt, InEop[i].yrt, eop.yrt - InEop[i].yrt      );
    */
    assert(eop.mjd.mjd() == InEop[i].mjd);
    assert(eop.xp == InEop[i].xp);
    assert(eop.yp == InEop[i].yp);
    assert(eop.dut == InEop[i].dut1);
    assert(eop.lod == InEop[i].lod);
    assert(eop.dx == InEop[i].dx);
    assert(eop.dy == InEop[i].dy);
    assert(eop.xrt == InEop[i].xrt);
    assert(eop.yrt == InEop[i].yrt);
  }

  return 0;
}
