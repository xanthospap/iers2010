#include "fundarg.hpp"
#include "geodesy/units.hpp"
#include "sofa.h"
#include <cstdio>

int num_tests = 10'000;
// constexpr const double MAX_ARCSEC = 1e-6;

int main() {

  double MAX_ARCSEC = 1e-6;
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random(40587, 66154);
    /* my result */
    const double mine = dso::fal03(mjd);
    /* SOFA result */
    const double sofa = iauFal03(mjd.jcenturies_sinceJ2000());
    /* assert difference < 1e-9 arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "fal03", MAX_ARCSEC, "OK");
 
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::falp03(mjd);
    /* SOFA result */
    const double sofa = iauFalp03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "falp03", MAX_ARCSEC, "OK");
  
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::faf03(mjd);
    /* SOFA result */
    const double sofa = iauFaf03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "faf03", MAX_ARCSEC, "OK");
  
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::fad03(mjd);
    /* SOFA result */
    const double sofa = iauFad03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    printf("[SOFA] %+.15f\n", sofa);
    printf("       %+.15f\n", mine);
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "fad03", MAX_ARCSEC, "OK");

  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::faom03(mjd);
    /* SOFA result */
    const double sofa = iauFaom03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "faom03", MAX_ARCSEC, "OK");
  
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::fame03(mjd);
    /* SOFA result */
    const double sofa = iauFame03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "fame03", MAX_ARCSEC, "OK");
  
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::fave03(mjd);
    /* SOFA result */
    const double sofa = iauFave03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "fave03", MAX_ARCSEC, "OK");
  
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::fae03(mjd);
    /* SOFA result */
    const double sofa = iauFae03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "fae03", MAX_ARCSEC, "OK");
  
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::fama03(mjd);
    /* SOFA result */
    const double sofa = iauFama03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "fama03", MAX_ARCSEC, "OK");
  
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::faju03(mjd);
    /* SOFA result */
    const double sofa = iauFaju03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "faju03", MAX_ARCSEC, "OK");
  
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::fasa03(mjd);
    /* SOFA result */
    const double sofa = iauFasa03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "fasa03", MAX_ARCSEC, "OK");
  
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::faur03(mjd);
    /* SOFA result */
    const double sofa = iauFaur03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "faur03", MAX_ARCSEC, "OK");
  
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::fane03(mjd);
    /* SOFA result */
    const double sofa = iauFane03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "fane03", MAX_ARCSEC, "OK");
  
  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random();
    /* my result */
    const double mine = dso::fapa03(mjd);
    /* SOFA result */
    const double sofa = iauFapa03(mjd.jcenturies_sinceJ2000());
    /* assert difference in arcsec */
    assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
  printf("%.20s %.3e %.5s\n", "fapa03", MAX_ARCSEC, "OK");
  
  return 0;
}
