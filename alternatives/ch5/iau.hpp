#ifndef __IAU_IERS10_SOFA_CPP_HPP__
#define __IAU_IERS10_SOFA_CPP_HPP__

#include "iersc.hpp"
#include <cmath>

namespace iers2010::sofa {

/// @brief Fundamental argument, IERS Conventions (2003): mean anomaly of the 
///        Moon.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT, 
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is 
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean anomaly of the Moon in radians
double fal03(double t) noexcept {
  double a =
      std::fmod(485868.249036 +
                    t * (1717915923.2178 +
                         t * (31.8792 + t * (0.051635 + t * (-0.00024470)))),
                TURNAS) *
      DAS2R;
  return a;
}

/// @brief Fundamental argument, IERS Conventions (2003): mean anomaly of the 
///        Sun.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT, 
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is 
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean anomaly of the Sun in radians
double falp03(double t) noexcept {
  double a =
      std::fmod(1287104.793048 +
                    t * (129596581.0481 +
                         t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))),
                TURNAS) *
      DAS2R;
  return a;
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of the 
///        Moon minus mean longitude of the ascending node
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT, 
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is 
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean anomaly of the Sun in radians
double fafp03(double t) noexcept {
  double a =
      std::fmod(335779.526232 +
                    t * (1739527262.8478 +
                         t * (-12.7512 + t * (-0.001037 + t * (0.00000417)))),
                TURNAS) *
      DAS2R;
  return a;
}

/// @brief Fundamental argument, IERS Conventions (2003): mean elongation of 
///        the Moon from the Sun.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT, 
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is 
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean anomaly of the Sun in radians
double fadp03(double t) noexcept {
  double a = std::fmod(          1072260.703692 +
             t * ( 1602961601.2090 +
             t * (        - 6.3706 +
             t * (          0.006593 +
             t * (        - 0.00003169 ) ) ) ), TURNAS ) * DAS2R;
  return a;
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of the 
///        Moon's ascending node.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT, 
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is 
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return Omega, radians
double faom03(double t) noexcept {
double a = std::fmod(          450160.398036 +
             t * ( - 6962890.5431 +
             t * (         7.4722 +
             t * (         0.007702 +
             t * (       - 0.00005939 ) ) ) ), TURNAS ) * DAS2R;
  return a;
}

/// @brief Fundamental argument, IERS Conventions (2003): mean longitude of 
///        Mercury.
/// @note  - Though t is strictly TDB, it is usually more convenient to use TT, 
///         which makes no significant difference.
///        - The expression used is as adopted in IERS Conventions (2003) and is 
///          from Simon et al. (1994).
///        - function is adopted from IAU SOFA, release 2021-05-12
/// @param[in] t TDB, Julian centuries since J2000.0
/// @return mean longitude of Mercury, radians
double fame03(double t) noexcept {
  return std::fmod(4.402608842 + 2608.7903141574 * t, D2PI);
}

}// ngpt::sofa

#endif