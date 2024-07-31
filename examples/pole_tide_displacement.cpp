#include "eop.hpp"
#include "geodesy/core/crdtype_warppers.hpp"
#include "pole_tide.hpp"

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s [opoleloadcoefcmcor.txt] [eop]\n", argv[0]);
    return 1;
  }

  /* put coordinates in a vector, in geodetic form */
  std::vector<SphericalCrd> sta;
  {
    /* SVAC at Ny-Alesund, Norway*/
    CartesianCrd svac;
    svac.mv << 1201300.042, 251874.432, 6238000.308;
    /* THULE at Thule Air Force base, Greenland */
    CartesianCrd thule;
    thule.mv << 538110.515, -1389031.364, 6180994.514;
    /* DIOB, Greece*/
    CartesianCrd diob;
    diob.mv << 4595212.468, 2039473.691, 3912617.891;
    /* KOUROU, French Guiana */
    CartesianCrd kourou;
    kourou.mv << 3855260.440, -5049735.535, 563056.590;
    /* BELGRANO at Antarctica */
    CartesianCrd belb;
    belb.mv << 1106046.627, -763739.010, -6214243.195;
    /* TERRE-ADELIE at Antarctica */
    CartesianCrd adhc;
    adhc.mv << -1940878.515, 1628473.041, -5833723.413;

    sta.emplace_back(cartesian2spherical(svac));
    sta.emplace_back(cartesian2spherical(thule));
    sta.emplace_back(cartesian2spherical(diob));
    sta.emplace_back(cartesian2spherical(kourou));
    sta.emplace_back(cartesian2spherical(belb));
    sta.emplace_back(cartesian2spherical(adhc));
  }

  const char *names[6] = { "svac", "thub", "diob", "krwb", "belb", "adhc"};

  /* a vector to hold interpolated coeffs */
  std::vector<OceanPoleTideDesaiCoeffs> coef;
  if (get_desai_ocp_deformation_coeffs(argv[1], sta, coef)) {
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
    const char **name = &names[0];
    for (const auto &site : sta) {
      /* (solid) pole tide displacement */
      auto dr1 = PoleTide::deformation(tt, eop.xp(), eop.yp(),
                                       SphericalCrdConstView(site));

      /* ocean pole tide displacement */
      auto dr2 = OceanPoleTide::deformation(tt, eop.xp(), eop.yp(),
                                            SphericalCrdConstView(site), *cf);

      /**/
      printf("%s %.9f %.9e %.9e %.9e %.9e %.9e %.9e\n", *name++,
             tt.imjd() + tt.fractional_days(), dr1.x(), dr1.y(), dr1.z(),
             dr2.x(), dr2.y(), dr2.z());

      ++cf;
    }

    /* next date */
    tt.add_seconds(FractionalSeconds(30));
  }
}
