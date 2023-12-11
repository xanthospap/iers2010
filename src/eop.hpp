/** @file
 * Utilities and parses to handle Earth Orientation Parameters (EOP).
 */
#ifndef __DSO_IERS_PRODUCTS_EOP_HPP__
#define __DSO_IERS_PRODUCTS_EOP_HPP__

#include "datetime/calendar.hpp"
#include <vector>
#include <array>
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
  using vit = std::vector<EopRecord>::iterator;
  using cvit = std::vector<EopRecord>::const_iterator;
  static constexpr int MAX_POLY_INTERPOLATION_DEGREE = 9;
  std::vector<EopRecord> mvec{};
  cvit last_index{mvec.end()};
  mutable std::array<double, (MAX_POLY_INTERPOLATION_DEGREE+1)*2> work;//[(MAX_POLY_INTERPOLATION_DEGREE+1)*2];
  //double *scr1() noexcept {return work;}
  //double *scr2() noexcept {return work + MAX_POLY_INTERPOLATION_DEGREE+1;}

public:
  enum class EopInterpoationResult : char {
    OutOfBoundsPrior,
    OutOfBoundsLater,
    Linear,
    PolyDegreeDescreased,
    PolyDegreeRequested
  }; /* EopInterpoationResult */

  //EopSeries()
  //    : mvec{}, last_index{mvec.end()},
  //      work(new double[(MAX_POLY_INTERPOLATION_DEGREE + 1) * 2]){};
  //~EopSeries() {
  //  if (work)
  //    delete[] work;
  //}

  /** Clear all entries in the series */
  void clear() noexcept {return mvec.clear();}

  /** Get number of epochs/entries in Series */
  int num_entries() const noexcept {return mvec.size();}

  /** Reserve capacity (not actual size!) 
   *
   * @param[in] capacity Number of EopRecord to reserve capacity for.
   */
  void reserve(int capacity) noexcept {mvec.reserve(capacity);}

  /** Push back an EopRecord entry.
   *
   * Warning! The vector of EopRecord should always be chronologically ordered. 
   * If yous use this method be sure that the entry you are pushing is indeed 
   * later that then last current entry in the table.
   */
  void push_back(const EopRecord& r) noexcept {return mvec.push_back(r);}

  EopInterpoationResult interpolate(const MjdEpoch &t, EopRecord &eop,
                                    int order = 5) const noexcept;

  /** @brief Return an iterator to the first Eop record in mvec, such that 
   * t < record.t
   *
   * Note that this means that:
   * 1. If the first entry in the vector is returned, t is out-of-bounds, 
   *    i.e. prior to the first EOP record.
   * 2. If the one-past-the-end iterator is returned (i.e. mvec.end()), then 
   *    t is out-of-bounds, i.e. at a latter epoch than the last EOP entry.
   */
  vit upper_bound(const MjdEpoch &t) noexcept {
    return std::upper_bound(mvec.begin(), mvec.end(),
                            [&](const EopRecord &r) { return t < r.t(); });
  }
  cvit upper_bound(const MjdEpoch &t) const noexcept {
    return std::upper_bound(mvec.begin(), mvec.end(),
                            [&](const EopRecord &r) { return t < r.t(); });
  }
}; /* EopSeries */

namespace details {
/* An enumeration type to define EOP format/series */
enum class IersEopFormat {C0414, C0420};

/* @brief Read header and decide the format of the passed in C04 file
 *
 * @param[in] c04fn An IERS EOP file in either C04/14 or C04/20 format
 * @param[out] type Depending on the header, the type of the input file
 * @return Anything other than zero means that the file was not identified! In 
 *         this case, type should not be used.
 */
int choose_c04_series(const char *c04fn, IersEopFormat &type) noexcept;

/* @brief Parse EOP data from an C04/14 IERS eop file for a given time range.
 *
 * Given a C04/14 IERS eop file, parse it and store the EOP data for the range
 * [start_mjd, end_mjd) to the passed-in EopSeries instance.
 *
 * Note that the EopLookUptable instance will be cleared before we start the
 * parsing (all previous data stored by the instance, if any, will be lost).
 * The data are stored in the same units as in the input EOP file, that is
 * arcsec for xp, yp, dX and Dy and seconds for dut1 and lod. Dates are
 * stored as TT MJDs.
 *
 * Note that for C04/14 files, the EOP data for xrt and yrt are not available.
 * This function will set the corresponding elements to 0.
 *
 * See https://hpiers.obspm.fr/iers/eop/eopc04_14/eopc04_IAU2000.62-now and
 * https://hpiers.obspm.fr/iers/eop/eopc04_14/updateC04.txt
 *
 * @param[in] c04fn An IERS C04/14 EOP file
 * @param[in] start_mjd Starting MJD to collect data for, in TT (inclusive)
 * @param[in] end_mjd Last MJD of the interval to collect data for (exclusive,
 *                    aka data will not be collected for this date) in TT
 * @param[in] eoptable A EopSeries where the parsed data will be stored
 *            in. The instance will be cleared, hence any data stored in it
 *            before the call, will be lost.
 * @return At success 0 is returned. Else, something went wrong
 */
int parse_iers_C0414(const char *c04fn, const MjdEpoch &start_tt,
                     const MjdEpoch &end_tt, EopSeries &eops) noexcept;

/* Given a C04/20 IERS eop file, parse it and store the EOP data for the 
 * dates [start_mjd, end_mjd) to the passed-in EopSeries instance.
 *
 * Note that the EopSeries instance will be cleared before we start the
 * parsing (all previous data stored by the instance, if any, will be lost).
 * The data are stored in the same units as in the input EOP file, that is
 * arcsec for xp, yp, dX and Dy and seconds for dut1 and lod. Dates are
 * stored as TT MJDs.
 *
 * See also https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.txt and
 * https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now
 * 
 * @param[in] c04fn An IERS C04/20 EOP file
 * @param[in] start_mjd Starting MJD to collect data for, in TT (inclusive)
 * @param[in] end_mjd Last MJD of the interval to collect data for (exclusive,
 *                    aka data will not be collected for this date) in TT
 * @param[in] eoptable A EopSeries where the parsed data will be stored
 *            in. The instance will be cleared, hence any data stored in it
 *            before the call, will be lost.  
 * @return At success 0 is returned. Else, something went wrong
 */
int parse_iers_C0420(const char *c04fn, const MjdEpoch &start_tt,
                     const MjdEpoch &end_tt, EopSeries &eops) noexcept;
} /* namespace details */

inline int parse_iers_C04(const char *c04fn, const MjdEpoch &start_tt,
                          const MjdEpoch &end_tt, EopSeries &eops) noexcept {
  details::IersEopFormat type;

  if (details::choose_c04_series(c04fn, type))
    return 1;

  if (type == details::IersEopFormat::C0420)
    return details::parse_iers_C0420(c04fn, start_tt, end_tt, eops);
  else
    return details::parse_iers_C0414(c04fn, start_tt, end_tt, eops);
}
} /* namespace dso */

#endif
