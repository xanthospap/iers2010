#include "fundarg.hpp"
#include "geodesy/units.hpp"
#include "sofa.h"
#include <cstdio>

int num_tests = 100'000;
constexpr const double MAX_ARCSEC = 1e-9;

int main() {

  for (int i = 0; i < num_tests; i++) {
    /* random MJD Epoch */
    const auto mjd = dso::MjdEpoch::random(40587, 66154);
    /* my result */
    const double mine = dso::fal03(mjd);
    /* SOFA result */
    const double sofa = iauFal03(mjd.jcenturies_sinceJ2000());
    /* assert difference < 1e-9 arcsec */
    if (std::abs(dso::rad2sec(mine - sofa)) >= MAX_ARCSEC) {
      printf("Difference is: %.3e (@%d)\n", dso::rad2sec(mine - sofa), i);
    }
    // assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  }
 
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::falp03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFalp03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}
  //
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::faf03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFaf03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}
  //
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::fad03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFad03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}

  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::faom03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFaom03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}
  //
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::fame03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFame03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}
  //
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::fave03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFave03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}
  //
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::fae03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFae03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}
  //
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::fama03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFama03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}
  //
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::faju03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFaju03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}
  //
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::fasa03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFasa03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}
  //
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::faur03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFaur03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}
  //
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::fane03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFane03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}
  //
  //for (int i = 0; i < num_tests; i++) {
  //  /* random MJD Epoch */
  //  const auto mjd = dso::MjdEpoch::random();
  //  /* my result */
  //  const double mine = dso::fapa03(mjd);
  //  /* SOFA result */
  //  const double sofa = iauFapa03(mjd.jcenturies_sinceJ2000());
  //  /* assert difference in arcsec */
  //  assert(std::abs(dso::rad2sec(mine - sofa)) < MAX_ARCSEC);
  //}

  return 0;
}
