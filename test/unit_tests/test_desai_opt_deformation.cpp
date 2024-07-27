#include "pole_tide.hpp"
#include <cstdio>
#include <cstring>
#include <vector>
#include "geodesy/transformations.hpp"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [opoleloadcoefcmcor.txt]\n", argv[0]);
    return 1;
  }

  /* SVAC at Ny-Alesund, Norway*/
  dso::CartesianCrd svac;
  svac.mv << 1201300.042, 251874.432, 6238000.308;

  /* THULE at Thule Air Force base, Greenland */
  dso::CartesianCrd thule;
  thule.mv << 538110.515, -1389031.364, 6180994.514;

  /* DIOB, Greece*/
  dso::CartesianCrd diob;
  diob.mv << 4595212.468, 2039473.691, 3912617.891;

  /* KOUROU, French Guiana */
  dso::CartesianCrd kourou;
  kourou.mv << 3855260.440, -5049735.535, 563056.590;

  /* BELGRANO at Antarctica */
  dso::CartesianCrd belb;
  belb.mv << 1106046.627, -763739.010, -6214243.195;

  /* TERRE-ADELIE at Antarctica */
  dso::CartesianCrd adhc;
  adhc.mv << -1940878.515, 1628473.041, -5833723.413;

  /* put coordinates in a vector, in geodetic form */
  std::vector<dso::GeodeticCrd> sta;
  sta.emplace_back(dso::cartesian2geodetic<dso::ellipsoid::grs80>(svac));
  sta.emplace_back(dso::cartesian2geodetic<dso::ellipsoid::grs80>(svac));
  sta.emplace_back(dso::cartesian2geodetic<dso::ellipsoid::grs80>(thule));
  sta.emplace_back(dso::cartesian2geodetic<dso::ellipsoid::grs80>(diob));
  sta.emplace_back(dso::cartesian2geodetic<dso::ellipsoid::grs80>(svac));
  sta.emplace_back(dso::cartesian2geodetic<dso::ellipsoid::grs80>(kourou));
  sta.emplace_back(dso::cartesian2geodetic<dso::ellipsoid::grs80>(belb));
  sta.emplace_back(dso::cartesian2geodetic<dso::ellipsoid::grs80>(adhc));
  sta.emplace_back(dso::cartesian2geodetic<dso::ellipsoid::grs80>(svac));

  /* a vector to hold interpolated coeffs */
  std::vector<dso::OceanPoleTideDesaiCoeffs> coef;

  if (dso::get_desai_ocp_deformation_coeffs(argv[1], sta, coef)) {
    fprintf(stderr, "ERROR Failed interpolating coeffs\n");
    return 1;
  }

  return 0;
}
