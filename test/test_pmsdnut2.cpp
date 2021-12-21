#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

constexpr double _alg_accuracy_ = 1e-6;

int main() {

  // testing pmsdnut2
  std::cout << "----------------------------------------\n";
  std::cout << "> pmsdnut2\n";
  std::cout << "----------------------------------------\n";

  const double dx_ref = 24.83144238273364834e0;
  const double dy_ref = -14.09240692041837661e0;
  const double fortran_diffs[] = {0.14e-06, 0.45e-06};

  double dx, dy;
  iers2010::pmsdnut2(54335e0, dx, dy);

  printf("\ndx= %12.6e microarcseconds", std::abs(dx - dx_ref));
  printf("\ndy= %12.6e microarcseconds", std::abs(dy - dy_ref));
  assert(std::abs(dx - dx_ref) < _alg_accuracy_ &&
         std::abs(dy - dy_ref) < _alg_accuracy_);
  assert(std::abs(dx - dx_ref) < fortran_diffs[0] && std::abs(dy - dy_ref) < fortran_diffs[1]);

  // check also the implementation using a datetime instance
  double dx2, dy2;
  dso::datetime<dso::seconds> t{dso::year(2007), dso::month(8),
                                dso::day_of_month(23), dso::seconds(0)};
  iers2010::pmsdnut2(t, dx2, dy2);
  assert(dx == dx2 && dy == dy2);

  dso::datetime<dso::seconds> t2(dso::modified_julian_day(54335e0),
                                 dso::seconds(0));
  iers2010::pmsdnut2(t2, dx2, dy2);
  assert(dx == dx2 && dy == dy2);

  return 0;
}
