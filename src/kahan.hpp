#ifndef __DSO_KAHAN_SUMMATION_ALG_HPP__
#define __DSO_KAHAN_SUMMATION_ALG_HPP__ 

#include <cmath>

namespace dso {

class KahanSum {
private:
  /** Accumulator */
  double msum{0e0};
  /** Running compensation for lost low-order bits */
  double merr{0e0};
public:
  /* improved Kahan–Babuška algorithm, see
   * https://en.wikipedia.org/wiki/Kahan_summation_algorithm
   */
  KahanSum &operator+=(double v) noexcept {
    const double t = msum + v;
    if (std::abs(msum) >= std::abs(v)) {
      /* If sum is bigger, low-order digits of v are lost. */
      merr += (msum-t) + v;
    } else {
      /* Else low-order digits of sum are lost. */
      merr += (v-t) + msum;
    }
    msum = t;
    /* return instance */
    return *this;
  }

  /** cast to double (after summation) */
  explicit operator double() const noexcept {return msum + merr;}

}; /* KahanSum */

} /* namespace dso */
#endif
