#include "aod1b.hpp"
#include "datetime/datetime_write.hpp"
#include <cstdio>

constexpr const dso::nanoseconds _01_hour_nanosec (60L     * 60L * 1'000'000'000L);
constexpr const dso::nanoseconds _03_hour_nanosec (3 * 60L * 60L * 1'000'000'000L);
constexpr const dso::nanoseconds _07_hour_nanosec (7 * 60L * 60L * 1'000'000'000L);
constexpr const dso::nanoseconds _23_hour_nanosec (23 * 60L* 60L * 1'000'000'000L);

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [AOD1B]\n", argv[0]);
    return 1;
  }

  dso::Aod1bIn aod(argv[1]);

  /* print some of the header stuff */
  printf("%.14e %.14e %.14e %.14e ", aod.GM(), aod.Re(), aod.flat(), aod.omega());

  /* store coefficients here */
  dso::StokesCoeffs cs(aod.max_degree(), aod.max_degree(), 0e0, 0e0);

  /* iterator for the ATM coefficients */
  dso::Aod1bDataBlockIterator<dso::AOD1BCoefficientType::ATM> it(aod);

  /* set iterator to the first data block */
  it.set_begin();

  /* benchmark epochs */
  auto bt1 = aod.first_epoch();
  bt1.add_seconds(_01_hour_nanosec);
  auto bt2 = aod.first_epoch();
  bt2.add_seconds(_03_hour_nanosec);
  auto bt3 = aod.first_epoch();
  bt3.add_seconds(_07_hour_nanosec);
  auto bt41 = aod.first_epoch();
  bt41.add_seconds(_23_hour_nanosec - _01_hour_nanosec - _01_hour_nanosec);
  auto bt42 = aod.first_epoch();
  bt42.add_seconds(_23_hour_nanosec);

  int record=1;
  int j = 0;
  while (!j) {
    if (it.header().mepoch <= bt1) {
      /* this should print coeffs of the first block */
      it.collect(cs);
      printf("%.14e %.14e ", cs.C(0,0), cs.S(0,0));
      printf("%.14e %.14e ", cs.C(10,0), cs.S(10,0));
      printf("%.14e %.14e ", cs.C(180,180), cs.S(180,180));
    } else if (it.header().mepoch <= bt2) {
      /* this should print coeffs of the second block */
      it.collect(cs, 120, 100);
      printf("%.14e %.14e ", cs.C(0,0), cs.S(0,0));
      printf("%.14e %.14e ", cs.C(10,0), cs.S(10,0));
      printf("%.14e %.14e ", cs.C(120,100), cs.S(120,100));
      printf("%.14e %.14e ", cs.C(120,101), cs.S(120,101));
    } else if (it.header().mepoch <= bt3) {
      /* this should print coeffs of the third block */
      it.collect(cs, 100, 100);
      printf("%.14e %.14e ", cs.C(0,0), cs.S(0,0));
      printf("%.14e %.14e ", cs.C(10,0), cs.S(10,0));
      printf("%.14e %.14e ", cs.C(100,100), cs.S(100,100));
      printf("%.14e %.14e ", cs.C(120,100), cs.S(120,100));
    } else if (it.header().mepoch >= bt41 && it.header().mepoch < bt42) {
      /* this should print coeffs of the last block */
      it.collect(cs, 180, 180);
      printf("%.14e %.14e ", cs.C(0,0), cs.S(0,0));
      printf("%.14e %.14e ", cs.C(10,0), cs.S(10,0));
      printf("%.14e %.14e ", cs.C(100,100), cs.S(100,100));
      printf("%.14e %.14e ", cs.C(120,100), cs.S(120,100));
    } else {
      it.skip();
    }
    ++record;
    j = it.advance();
  }

  /* should have reached EOF */
  assert(j<0);
  //printf("\n");

  return 0;
}
