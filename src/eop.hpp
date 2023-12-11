/** @file
 * Utilities and parses to handle Earth Orientation Parameters (EOP).
 */
#ifndef __DSO_IERS_PRODUCTS_EOP_HPP__
#define __DSO_IERS_PRODUCTS_EOP_HPP__

#include "datetime/calendar.hpp"
#include <vector>
#include <algorithm>

namespace dso {

/** @brief A structure to hold EOP records for a single epoch */
class EopRecord {
  dso::MjdEpoch mt;
  /* Pole coordinates x, y (") */
  double mxp, myp;
  /* UT1 - UTC [sec] */
  double mdut;
  /* Lenght of day offset LOD [sec] */
  double mlod;
  /* Celestial pole offsets dX, dY ["]; see last paragraphs of IERS 2010, 
   * Section 5.5.4 
   */
  double mdx, mdy;
  /* ΔAT at time mt */
  double mdat;
  /* Pole coordinate rates ["/day], only available in C04/20 */
  double mxrt, myrt;

public:
  /** Date in TODO */
  dso::MjdEpoch t() const noexcept {return mt;}
  dso::MjdEpoch &t() noexcept {return mt;}

  /* ΔAT at time mt */
  double dat() const noexcept {return mdat;}
  double &dat() noexcept {return mdat;}
  /** Pole coordinate xp in [sec] */
  double xp() const noexcept {return mxp;}
  /** Pole coordinate yp in [sec] */
  double yp() const noexcept {return myp;}
  /** UT1 - UTC in [sec] */
  double dut() const noexcept {return mdut;}
  /** Lenght of day offset LOD [sec] */
  double lod() const noexcept {return mlod;}
  /** Celestial pole offsets δX in [sec] */
  double dX() const noexcept {return mdx;}
  /** Celestial pole offsets δY in [sec] */
  double dY() const noexcept {return mdy;}
  /** Pole rate in X, in [sec/day] */
  double xp_rate() const noexcept {return mxrt;}
  /** Pole rate in Y, in [sec/day] */
  double yp_rate() const noexcept {return myrt;}
  /** Pole coordinate xp in [sec] */
  double &xp() noexcept {return mxp;}
  /** Pole coordinate yp in [sec] */
  double &yp() noexcept {return myp;}
  /** UT1 - UTC in [sec] */
  double &dut() noexcept {return mdut;}
  /** Lenght of day offset LOD [sec] */
  double &lod() noexcept {return mlod;}
  /** Celestial pole offsets δX in [sec] */
  double &dX() noexcept {return mdx;}
  /** Celestial pole offsets δY in [sec] */
  double &dY() noexcept {return mdy;}
  /** Pole rate in X, in [sec/day] */
  double &xp_rate() noexcept {return mxrt;}
  /** Pole rate in Y, in [sec/day] */
  double &yp_rate() noexcept {return myrt;}

  /* @brief Angular velocity of Earth in [rad/sec], including LOD variation
   * @return Angular velocity of Earth (Ω) in [rad/sec]
   */
  //double omega() const noexcept {
  //  /* see https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html */
  //  const double OmPicoRadSec =
  //      iers2010::OmegaEarth * 1e12 - 0.843994809e0 * (lod * 1e3);
  //  return OmPicoRadSec * 1e-12;
  //}
};/* EopRecord */

class EopSeries {
  std::vector<EopRecord> mvec;
  using vit = std::vector<EopRecord>::iterator;
  using cvit = std::vector<EopRecord>::const_iterator;

public:
  enum class EopInterpoationResult : char {
    OutOfBoundsPrior,
    OutOfBoundsLater,
    Linear,
    PolyDegreeDescreased,
    PolyDegreeRequested
  }; /* EopInterpoationResult */

  EopInterpoationResult interpolate(const MjdEpoch &t, EopRecord &eop,
                                    int order = 5) const;

  /** @brief Return an iterator to the first Eop record in mvec, such that 
   * t < record.t
   *
   * Note that this means that:
   * 1. If the first entry in the vector is returned, t is out-of-bounds, 
   *    i.e. prior to the first EOP record.
   * 2. If the one-past-the-end iterator is returned (i.e. mvec.end()), then 
   *    t is out-of-bounds, i.e. at a latter epoch than the last EOP entry.
   */
  vit upper_bound(const MjdEpoch &t) const noexcept {
    return std::upper_bound(mvec.begin(), mvec.end(),
                            [&](const EopRecord &r) { return t < r.t(); });
  }
}; /* EopSeries */
} /* namespace dso */

#endif
