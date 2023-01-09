#include "doodson.hpp"
#include "hardisp.hpp"
#include "iers2010.hpp"
#include <datetime/dtcalendar.hpp>
#ifdef DEBUG
#include <cassert>
#endif

namespace {
Constituent Constituents[] = {
    {{2, 5, 5, 5, 5, 5}, +6.3221e-01},  {{2, 7, 3, 5, 5, 5}, +2.9411e-01},
};
constexpr const int NT = sizeof(Constituents) / sizeof(Constituents[0]);

}// unnamed namespace

class Hardisp {
  double beta[6], beta_freq[6];

  int operator()(const dso::TwoPartDate &tt_mjd,
                 const dso::TwoPartDate &ut1_mjd) noexcept {
    const dso::TwoPartDate tt_jd = tt_mjd.jd_sofa();
    const dso::TwoPartDate ut1_jd = ut1_mjd.jd_sofa();
    // compute angles:
    // 1. Fundamental Arguments (temporary)
    double fundarg[5];
    iers2010::fundarg(tt_mjd.jcenturies_sinceJ2000(), fundarg);
    // 2. GMST [rad]
    const double gmst = iers2010::sofa::gmst06(ut1_jd._big, ut1_jd._small,
                                               tt_jd._big, tt_jd._small);
    // 3. Doodson variables (stored in beta)
    dso::fundarg2doodson(fundarg, gmst, beta);
    // 4. Variables for frequency (storedin in beta_freq)
    dso::DoodsonNumber::doodson_freq_vars(tt_mjd, beta_freq);

    return 0;
  }
};

int foo(const dso::TwoPartDate &tt_mjd, const dso::TwoPartDate &ut1_mjd,
        const std::vector<Constituent> &cin) noexcept {
  // compute angles
  this->operator()(tt_mjd, ut1_mjd);

  for (const auto &i : cin) {
    // match given constituent
    auto it = std::find_if(Constituents.begin(), Constituents.end(),
                           [&i = std::as_const(i)](const Constituent &c) {
                             return c.doodson == j.doodson;
                           });
    if (it != Constituents.end()) {
      const double phase = i.phase(beta);
      const double freq = i.frequency(beta_freq);
    }
  }
}
