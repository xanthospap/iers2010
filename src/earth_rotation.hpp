/** @file
 * Class and unctions to implement an easy-to-use API to describe Earth
 * rotation and the ties between ITRS and GCRS frames.
 */

#ifndef __DSO_EARTH_ROTATION_API_HPP__
#define __DSO_EARTH_ROTATION_API_HPP__

#include "eop.hpp"
#include "eigen3/Eigen/Eigen"
#include "eigen3/Eigen/src/Geometry/Quaternion.h"

namespace dso {
class EarthRotation {
private:
  /** EOPs to use for interpolation/trnaformation */
  EopSeries meops;
  /** Last used EOP values */
  EopRecord meop;

  int eops_at(const MjdEpoch &tt) noexcept;

public:
  /** Constructor using a C04 EOP file (IERS) */
  EarthRotation(const char *fn, const MjdEpoch &start_tt = MjdEpoch::min(),
                const MjdEpoch &end_tt = MjdEpoch::max());

  /** Get the unit quaternion to transform from ITRS to GCRS */
  Eigen::Quaterniond itrs2gcrs_quaternion(const MjdEpoch &tt,
                                          double *fargs14=nullptr);
}; /* class EarthRotation */
} /* namespace dso */

#endif
