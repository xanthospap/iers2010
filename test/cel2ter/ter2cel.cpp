#include "cel2ter.hpp"
#include <datetime/dtcalendar.hpp>

// Jason 3 orbit in ITRF
struct Sp3Orb {
  dso::datetime<dso::nanoseconds> t_;
  Eigen::Matrix<double,3,1> pos_,vel_;
  dso::TwoPartDate tai() const noexcept {return dso::TwoPartDate(t_);}
  Eigen::Matrix<double,3,1> pos() const noexcept {return pos_;}
  Eigen::Matrix<double,3,1> vel() const noexcept {return vel_;}
};

const Sp3Orb Data[] = {
    {dso::datetime<dso::nanoseconds>(dso::year(2020), dso::month(12),
                                     dso::day_of_month(27), dso::hours(8),
                                     dso::minutes(1), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({3487.293316, -896.788983, -6824.852203}),
     Eigen::Matrix<double, 3, 1>({-15363.302305, 65807.746967, -16487.845173})},
    {dso::datetime<dso::nanoseconds>(dso::year(2020), dso::month(12),
                                     dso::day_of_month(28), dso::hours(13),
                                     dso::minutes(14), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({-3515.619511, 1582.273460, 6682.500148}),
     Eigen::Matrix<double, 3, 1>(
         {-55881.071368, -35814.095603, -20898.293109})},
    {dso::datetime<dso::nanoseconds>(dso::year(2020), dso::month(12),
                                     dso::day_of_month(29), dso::hours(10),
                                     dso::minutes(1), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({-3742.495386, -5124.360763, 4390.586877}),
     Eigen::Matrix<double, 3, 1>({3932.260187, -46937.326227, -51379.107782})},
    {dso::datetime<dso::nanoseconds>(dso::year(2020), dso::month(12),
                                     dso::day_of_month(30), dso::hours(8),
                                     dso::minutes(45), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({209.489385, -7635.120039, -1108.019286}),
     Eigen::Matrix<double, 3, 1>({23643.552913, 10061.562733, -64878.143357})},
    {dso::datetime<dso::nanoseconds>(dso::year(2020), dso::month(12),
                                     dso::day_of_month(30), dso::hours(17),
                                     dso::minutes(34), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({2887.269309, 1498.511837, 6995.113390}),
     Eigen::Matrix<double, 3, 1>({-46306.375811, 51275.925785, 8123.030716})},
    {dso::datetime<dso::nanoseconds>(dso::year(2020), dso::month(12),
                                     dso::day_of_month(31), dso::hours(14),
                                     dso::minutes(51), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({-6713.956408, -789.405749, -3723.422442}),
     Eigen::Matrix<double, 3, 1>({33838.228849, -24592.150829, -55773.754675})},
    {dso::datetime<dso::nanoseconds>(dso::year(2021), dso::month(1),
                                     dso::day_of_month(1), dso::hours(23),
                                     dso::minutes(53), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({340.779362, -3430.909701, 6901.240056}),
     Eigen::Matrix<double, 3, 1>({59834.266448, 32880.895867, 13379.748152})},
    {dso::datetime<dso::nanoseconds>(dso::year(2021), dso::month(1),
                                     dso::day_of_month(2), dso::hours(23),
                                     dso::minutes(31), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({-4611.654958, 104.520344, -6186.375718}),
     Eigen::Matrix<double, 3, 1>({-43284.998727, -44479.014894, 31497.334336})},
    {dso::datetime<dso::nanoseconds>(dso::year(2021), dso::month(1),
                                     dso::day_of_month(4), dso::hours(5),
                                     dso::minutes(56), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({-4114.785663, 5951.943373, 2682.293534}),
     Eigen::Matrix<double, 3, 1>({-8435.969148, -33240.371176, 60748.068472})},
    {dso::datetime<dso::nanoseconds>(dso::year(2021), dso::month(1),
                                     dso::day_of_month(4), dso::hours(17),
                                     dso::minutes(53), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({-1900.968447, 7007.042604, 2616.092961}),
     Eigen::Matrix<double, 3, 1>({-30583.139278, 14502.925760, -60996.445094})},
    {dso::datetime<dso::nanoseconds>(dso::year(2021), dso::month(1),
                                     dso::day_of_month(5), dso::hours(2),
                                     dso::minutes(9), dso::nanoseconds(0L)),
     Eigen::Matrix<double, 3, 1>({-4033.082313, 3326.742862, -5676.125130}),
     Eigen::Matrix<double, 3, 1>({-57618.710859, -3358.060832, 38949.778141})}};

const int DataSize = sizeof(Data) / sizeof(Data[0]);

int main(int argc, char *argv[]) {
  if (argc !=2 ) {
    fprintf(stderr, "Usage: %s [C04/14 EOP file] \n", argv[0]);
    return 1;
  }

  // Create an EOP LookUp table
  dso::EopLookUpTable eop_lut;
  {
    const int start = (int)Data[0].tai().mjd() - 5;
    const int end = (int)Data[DataSize - 1].tai().mjd() + 5;
    // parse C04 EOPs
    if (dso::parse_iers_C0414(argv[1], start, end, eop_lut)) {
      fprintf(stderr, "ERROR. Failed collecting EOP data\n");
      return 1;
    }
    eop_lut.utc2tt();
    eop_lut.regularize(false);
  }

  // one by one, transform sp3 input orbit from ITRF to GCRF and back
  // check results
  for (int i=0; i<DataSize; i++) {
    dso::Itrs2Gcrs Rot(Data[i].tai().tai2tt(), &eop_lut);
    const Eigen::Matrix<double,3,1> r_gcrf = Rot.itrf2gcrf(Data[i].pos());
    const Eigen::Matrix<double,3,1> r_itrf = Rot.gcrf2itrf(r_gcrf);
    printf("ITRF Diffs: %+.9f %+.9f %+.9f\n", Data[i].pos()(0) - r_itrf(0),
           Data[i].pos()(1) - r_itrf(1), Data[i].pos()(2) - r_itrf(2));
  }

  return 0;
}
