#include "pole_tide.hpp"
#include <cstdio>
#include <cstring>
#include <vector>
#include "geodesy/transformations.hpp"
#include "geodesy/units.hpp"
#include <cassert>

std::vector<dso::OceanPoleTideDesaiCoeffs> reference_results = {
dso::OceanPoleTideDesaiCoeffs{-0.083842057482,-0.031487496025,-0.052192091856,-0.049714846909,+0.017543805543,+0.063145164562},
dso::OceanPoleTideDesaiCoeffs{-0.083842057482,-0.031487496025,-0.052192091856,-0.049714846909,+0.017543805543,+0.063145164562},
dso::OceanPoleTideDesaiCoeffs{-0.042708879444,+0.004020731856,+0.007115063343,+0.035590298032,+0.046970420255,+0.047096934154},
dso::OceanPoleTideDesaiCoeffs{-0.209475806381,-0.060027439515,-0.044124236270,-0.065587650505,+0.005072914900,+0.058382406965},
dso::OceanPoleTideDesaiCoeffs{-0.083842057482,-0.031487496025,-0.052192091856,-0.049714846909,+0.017543805543,+0.063145164562},
dso::OceanPoleTideDesaiCoeffs{-0.023511531496,+0.007228134473,-0.017765274488,-0.040658975747,+0.046938342434,+0.062856621739},
dso::OceanPoleTideDesaiCoeffs{+0.078336289037,+0.005541119869,-0.001750128778,-0.014520796339,+0.030240602665,+0.062190940126},
dso::OceanPoleTideDesaiCoeffs{-0.123011087288,+0.163345531723,-0.007347183834,-0.001903472227,-0.031746368694,-0.059410295632},
dso::OceanPoleTideDesaiCoeffs{-0.083842057482,-0.031487496025,-0.052192091856,-0.049714846909,+0.017543805543,+0.063145164562}};

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

  auto cr = reference_results.cbegin();
  for (const auto &c : coef) {
    // printf("%.9f %.9f %.9f %.9f %.9f %.9f\n", c.rR,c.rI,c.nR,c.nI,c.eR,c.eI);
    // printf("%.9f %.9f %.9f %.9f %.9f %.9f\n", cr->rR,cr->rI,cr->nR,cr->nI,cr->eR,cr->eI);
    assert(std::abs(c.rR - cr->rR) < 1e-9);
    assert(std::abs(c.rI - cr->rI) < 1e-9);
    assert(std::abs(c.nR - cr->nR) < 1e-9);
    assert(std::abs(c.nI - cr->nI) < 1e-9);
    assert(std::abs(c.eR - cr->eR) < 1e-9);
    assert(std::abs(c.eI - cr->eI) < 1e-9);
    ++cr;
  }

  return 0;
}
