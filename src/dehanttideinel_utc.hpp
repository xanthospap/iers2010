#ifndef __DEHANTTIDEINEL_UTC_IMPL_HPP__
#define __DEHANTTIDEINEL_UTC_IMPL_HPP__

#include "iers2010.hpp"

namespace iers2010 {

#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
Eigen::Matrix<double, 3, 1>
dehanttideinel_from_utc(const Eigen::Matrix<double, 3, 1>& xsta,
    const Eigen::Matrix<double, 3, 1>& xsun,
    const Eigen::Matrix<double, 3, 1>& xmon,
    const dso::datetime<S>& tutc) noexcept
{
  auto t = tutc;
  
  // need to convert UTC to TT: TT = UTC + 32.184[sec] + DAT = 
  int dat = dso::dat(t.mjd());
  t.add_seconds(dso::milliseconds(dat * 1e3) + dso::milliseconds(32184));
  
  return dehanttideinel_impl(xsta, xsun, xmon, t.jcenturies_sinceJ2000(),
      tutc.sec().to_fractional_seconds() / 3600e0);
}

} // iers2010
#endif
