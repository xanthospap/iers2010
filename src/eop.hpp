/** @file
 * Utilities and parsers to handle Earth Orientation Parameters (EOP) and 
 * (time)series of such values.
 */
#ifndef __DSO_IERS_PRODUCTS_EOP_HPP__
#define __DSO_IERS_PRODUCTS_EOP_HPP__

#include "datetime/calendar.hpp"
#include <vector>
#include <array>
#include <algorithm>
#include <utility>

namespace dso {

/** Compute libration effect, i.e. the diurnal lunisolar effect on polar 
 *  motion.
 *
 * The libration effect is described in IERS 2010, Ch. 5.5.1.3 "Variations 
 * (∆x, ∆y)libration in polar motion". 
 * This functions performs the computations published in the PMSDNUT2 
 * subroutine available via the IERS 
 * https://iers-conventions.obspm.fr/content/chapter5/software/PMSDNUT2.F
 *
 * @param[in] fargs The fundamental arguments
 * @param[in] gmst  The GMST
 * @param[out] xp   Δx due to libration in x-pole [microarcseconds]
 * @param[out] yp   Δy due to libration in y-pole [microarcseconds]
 * @return Always 0
 */
int xy_libration(const double *const fargs, double gmst, double &xp,
                 double &yp) noexcept;
int ut_libration(const double *const fargs, double gmst, double &dut1,
                     double &dlod) noexcept;

/** @brief A structure to hold EOP records for a single epoch */
class EopRecord {
  /** Record epoch in TT */
  dso::MjdEpoch mt;
  /** Pole coordinates x, y (") */
  double mxp, myp;
  /** UT1 - UTC [sec] */
  double mdut;
  /** Lenght of day offset LOD [sec] */
  double mlod;
  /* Celestial pole offsets dX, dY ["]; see last paragraphs of IERS 2010, 
   * Section 5.5.4 
   */
  double mdx, mdy;
  /** ΔAT at time mt */
  double mdat;
  /** Pole coordinate rates ["/day], only available in C04/20 */
  double mxrt, myrt;

public:
  /** Date in TT */
  dso::MjdEpoch t() const noexcept {return mt;}
  /** Date in TT */
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
  /** Maximum allowed interpolation degree */
  static constexpr int MAX_POLY_INTERPOLATION_DEGREE = 9;
  /** Vector of EopRecords stored in series */
  std::vector<EopRecord> mvec{};
  /** Last used index in mvec; speed-up interpolation when requesting 
   * interpolation results for epochs within one day (e.g. every N seconds). 
   * This is the most often use-case.
   * This is a mutable member, i.e. changing values in here does not actually 
   * 'change' the instance; it is only used to hunt interpolation indexes.
   */
  mutable cvit last_it{mvec.end()};
  /** Scratch/Workspace memoty to use in interpolation. This is a mutable 
   * member, i.e. changing values in here does not actually 'change' the 
   * instance.
   */
  mutable std::array<double, (MAX_POLY_INTERPOLATION_DEGREE+1)*2> work;

public:

  /** Signal interpolation result */
  enum class EopInterpolationResult : char {
    OutOfBoundsPrior,
    OutOfBoundsLater,
    Linear,
    PolyDegreeDescreased,
    PolyDegreeRequested
  }; /* EopInterpolationResult */

  const auto &give_me_the_vector() const noexcept {
    return mvec;
  }

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
   *
   * Warning! This function may invalidate the last_it pointer. Please, reset 
   * it to the correct index after the call.
   */
  void push_back(const EopRecord& r) noexcept {return mvec.push_back(r);}

  void reset_vec_iterator(int offset=-1) noexcept {
    last_it = (offset < 0) ? mvec.end() : mvec.begin()+offset;
  }

  /** Interpolate EOPs using the data in Series, for a given epoch (TT).
   *
   * a. In case \p t is prior to the first entry in the instant's series, i.e. 
   *    t <= mvec[0].t:
   *    \p eop will be filled with the values of the first entry in the series, 
   *    (i.e. no interpolation) aka eop == mvec[0] (apparently the equality 
   *    excludes the epoch of \p eop, which will be \t).
   *    The return status will be set to 
   *    EopSeries::EopInterpolationResult::OutOfBoundsPrior
   *
   * b. In case \p t is later than the last entry in the instant's series, i.e. 
   *    t >= mvec[num_entries].t:
   *    \p eop will be filled with the values of the last entry in the series, 
   *    (i.e. no interpolation) aka eop == mvec[num_entries] (apparently the 
   *    equality excludes the epoch of \p eop, which will be \t).
   *    The return status will be set to 
   *    EopSeries::EopInterpolationResult::OutOfBoundsLater
   * 
   * c. For any other case, the function will try to use (degree+1)/2 data 
   *    points on the left and right-hand side (i.e. in total degree+1 data 
   *    points) to perform polynomial interpolation via the Lagrange method.
   *    If this is not possible, the function will test lowering the degree of 
   *    the interpolating polynomial (i.e. degree, degree-1, ...,1), untill it 
   *    can find a suitable range of (degree+1)/2 data points on the left and 
   *    right-hand side; it will then perform the interpolation using this 
   *    degree/num of data points. If the 'final' degree is lower than the 
   *    user-specified (i.e. \p degree), then the status 
   *    EopSeries::EopInterpolationResult::PolyDegreeDescreased will be 
   *    returned. In the special case where only one point on the left and one 
   *    on the right can be used (i.e. degree=1), the status function will 
   *    perform linear interpolation and return 
   *    EopSeries::EopInterpolationResult::Linear. If no decreasing of the 
   *    specified degree was needed, the function will return the status
   *    EopSeries::EopInterpolationResult::PolyDegreeRequested.
   *
   * Note that the function will always use the same number of data points on 
   * the left and right-hand side to perform interpolation. E.g., in case 
   * degree=5, three data points prior to \p t and three later than \p t.
   *
   * Interpolation for the ΔUT1 EOP parameter, is performed as descibed in:
   * Ben K. Bradley, Aurore Sibois, Penina Axelrad, Influence of ITRS/GCRS 
   * implementation for astrodynamics: Coordinate transformations, Advances in 
   * Space Research, Volume 57, Issue 3, 2016, Pages 850-866, ISSN 0273-1177,
   * https://doi.org/10.1016/j.asr.2015.11.006. That means that we:
   * 1. Construct continuous UT1-TAI data points,
   * 2. Interpolate at the requested epoch, and
   * 3. Add leap seconds for day in question
   * Hence, the value returned for ΔUT1 is "correctly" interpolated over leap 
   * seconds.
   *
   * @warning The ΔΑΤ parameter of the EOP records in the series, are excluded 
   *          from the interpolation; hence, the resulting EopRecord, will 
   *          have a random/invalid number of ΔΑΤ seconds. 
   *          **Do not use it**.
   *
   * @param[in] t Epoch at which the interpolation is requested, in TT.
   * @param[out] eop Interpolation result, for all EOPs, excluding ΔAT.
   * @param[in] degree Poynomial interpolation degree. The data points used 
   *            to interpolate will be degree+1; (degree+1)/2 on the left and 
   *            another (degree+1)/2 on the right. Maximum allowed degree, is
   *            MAX_POLY_INTERPOLATION_DEGREE.
   */
  EopInterpolationResult interpolate(const MjdEpoch &t, EopRecord &eop,
                                    int degree = 5) const noexcept;

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
    return std::upper_bound(
        mvec.begin(), mvec.end(), t,
        [/*&t = std::as_const(t)*/](const MjdEpoch &e, const EopRecord &r) {
          return e < r.t();
        });
  }
  cvit upper_bound(const MjdEpoch &t) const noexcept {
    return std::upper_bound(
        mvec.begin(), mvec.end(), t,
        [/*&t = std::as_const(t)*/](const MjdEpoch &e, const EopRecord &r) {
          return e < r.t();
        });
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

int parse_iers_C04(const char *c04fn, const MjdEpoch &start_tt,
                   const MjdEpoch &end_tt, EopSeries &eops) noexcept;
} /* namespace dso */

#endif
