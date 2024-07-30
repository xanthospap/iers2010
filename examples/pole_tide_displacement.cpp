#include "eop.hpp"
#include "pole_tide.hpp"
#include <geodesy/core/crdtype_warppers.hpp>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s [opoleloadcoefcmcor.txt] [eop]\n", argv[0]);
    return 1;
  }

  /* put coordinates in a vector, in geodetic form */
  std::vector<dso::SphericalCrd> sta;
  {
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

    sta.emplace_back(dso::cartesian2spherical(svac));
    sta.emplace_back(dso::cartesian2spherical(thule));
    sta.emplace_back(dso::cartesian2spherical(diob));
    sta.emplace_back(dso::cartesian2spherical(kourou));
    sta.emplace_back(dso::cartesian2spherical(belb));
    sta.emplace_back(dso::cartesian2spherical(adhc));
  }

  /* a vector to hold interpolated coeffs */
  std::vector<dso::OceanPoleTideDesaiCoeffs> coef;
  if (dso::get_desai_ocp_deformation_coeffs(argv[1], sta, coef)) {
    fprintf(stderr, "ERROR Failed interpolating coeffs\n");
    return 1;
  }

  /* start and end epochs in TT */
  MjdEpoch start(year(2023), month(7), day_of_month(19));
  MjdEpoch end(year(2023), month(7), day_of_month(29));

  /* create an instance to hold EOP series */
  EopSeries eops;

  /* parse EOP values to the EopSeries instance for any t >= t1 and t < t2 */
  if (parse_iers_C04(argv[2], start - modified_julian_day(2),
                     end + modified_julian_day(2), eops)) {
    fprintf(stderr, "ERROR Failed parsing eop file\n");
    return 1;
  }

  auto tt = start;
  while (tt < end) {
    /* EOPS at time of request */
    EopRecord eop;
    if (EopSeries::out_of_bounds(eops.interpolate(tt, eop))) {
      fprintf(stderr, "Failed to interpolate: Epoch is out of bounds!\n");
      return 1;
    }

    auto cf = coef.begin();
    for (const auto &site : sta) {
      /* (solid) pole tide displacement */
      auto dr1 = dso::PoleTide::deformation(tt, eop.xp(), eop.yp(), site);

      /* ocean pole tide displacement */
      auto dr2 =
          dso::OceanPoleTide::deformation(tt, eop.xp(), eop.yp(), site, cf);

      /**/
      printf("%.9f %.4f %.4f %.4f %.4f%.4f %.4f\n",
             tt.imjd() + tt.fractional_days(), dr1(0), dr1(1), dr1(2), dr2(0),
             dr2(1), dr2(2));

      ++cf;
    }

    /* next date */
    tt.add_seconds(FractionalSeconds(30));
  }
}
