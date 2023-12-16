/** @file
 * Define a class to handle the ITRS to GCRS (and the inverse) transformation, 
 * both for position and velocity vector
 */

#ifndef __DSO_ITRS_TO_GCRS_ROTATION_HPP__
#define __DSO_ITRS_TO_GCRS_ROTATION_HPP__

#include "eop.hpp"
#include "iau.hpp"

namespace dso {

class RotItrsGcrs {
private:
  /** An EopSeries instance, i.e. EOP data to compute rotation */
  EopSeries meopsv;
}; /* RotItrsGcrs */

} /* namespace dso */

#endif
