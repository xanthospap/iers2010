#include <cstring>
#include "doodson.hpp"
#include "iau.hpp"
#include "fundarg.hpp"

char *dso::DoodsonNumber::str(char *buf,
                              bool use_5s_convention) const noexcept {
  int offset = use_5s_convention ? 5 : 0;
  sprintf(buf, "%d%d%d.%d%d%d", iar[0], iar[1] + offset, iar[2] + offset,
          iar[3] + offset, iar[4] + offset, iar[5] + offset);
  return buf;
}

double dso::DoodsonNumber::phase(
    const dso::TwoPartDate &tt_mjd,
    const dso::TwoPartDate &ut1_mjd) const noexcept {
  const dso::TwoPartDate tt_jd = tt_mjd.jd_sofa();
  const dso::TwoPartDate ut1_jd = ut1_mjd.jd_sofa();
  // Greenwich mean sidereal time (radians)
  const double gmst = iers2010::sofa::gmst06(ut1_jd._big, ut1_jd._small,
                                             tt_jd._big, tt_jd._small);
  // Fundamental arguments
  double fargs[5], beta[6];
  iers2010::fundarg(tt_mjd.jcenturies_sinceJ2000(), fargs);
  // Doodson arguments (store in beta[0:6)
  dso::fundarg2doodson(fargs, gmst, beta);
  return phase(beta);
}

/// @brief Frequency of this constituent in [rad/day]
double *dso::DoodsonNumber::doodson_freq_vars(const dso::TwoPartDate &tt_mjd,
                                              double *args) noexcept {
  const double t = tt_mjd.jcenturies_sinceJ2000()/10e0;
  args[0] = dso::anp(dso::deg2rad(
      127037328.88553056e0 +
      (2 * 0.17696111e0 + (3 * -0.0018314e0 + 4 * 0.00008824e0 * t) * t) * t));
  args[1] = dso::anp(dso::deg2rad(
      4812678.8119575e0 +
      (2 * -0.14663889e0 + (3 * 0.0018514e0 + 4 * -0.00015355e0 * t) * t) * t));
  args[2] = dso::anp(dso::deg2rad(
      360007.69748806e0 +
      (2 * 0.03032222e0 + (3 * 0.00002e0 + 4 * -0.00006532e0 * t) * t) * t));
  args[3] = dso::anp(dso::deg2rad(
      40690.1363525e0 +
      (2 * -1.03217222e0 + (3 * -0.01249168e0 + 4 * 0.00052655e0 * t) * t) *
          t));
  args[4] = dso::anp(dso::deg2rad(
      19341.36261972e0 +
      (2 * -0.20756111e0 + (3 * -0.00213942e0 + 4 * 0.00016501e0 * t) * t) *
          t));
  args[5] = dso::anp(dso::deg2rad(
      17.19457667e0 +
      (2 * 0.04568889e0 + (3 * -0.00001776e0 + 4 * -0.00003323e0 * t) * t) *
          t));
  return args;
}
