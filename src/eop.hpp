#ifndef __DSO_IERS_PRODUCTS_EOP_HPP__
#define __DSO_IERS_PRODUCTS_EOP_HPP__

#include "iers2010.hpp"
#include <vector>
#ifdef DEBUG
#include <cassert>
#endif

namespace dso {

/* @bref A structure to hold EOP records for a single epoch */
struct EopRecord {
  dso::TwoPartDate mjd;
  /* Pole coordinates x, y (") */
  double xp, yp;
  /* UT1 - UTC [sec] */
  double dut;
  /* Lenght of day offset LOD [sec] */
  double lod;
  /* Celestial pole offsets dX, dY ["] */
  double dx, dy;
  /* Pole coordinate rates ["/day], only available in C04/20 */
  double xrt, yrt;

  /* @brief Angular velocity of Earth in [rad/sec], including LOD variation
   * @return Angular velocity of Earth (Ω) in [rad/sec]
   */
  double omega() const noexcept {
    /* see https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html */
    const double OmPicoRadSec =
        iers2010::OmegaEarth * 1e12 - 0.843994809e0 * (lod * 1e3);
    return OmPicoRadSec * 1e-12;
  }
};/* EopRecord */

/* @brief A class to hold EOP information ordered chronologically */
class EopLookUpTable {
  /* UTC or TT */
  std::vector<dso::TwoPartDate> tvec;
  std::vector<double> xp, yp, dut1, dX, dY, lod, xrt, yrt;

public:
  /* @brief Clear all instance's arrays */
  void clear() noexcept;

#ifdef DEBUG
  void printxp() const noexcept { 
    printf("Xp Series:\n");
    for (const auto x : xp) printf("%+9.6f, ", x);
    printf("\n");
  }
#endif
  
  /* @brief Return an iterator to the tvec vector, such that: it > t */
  auto upper_bound(const dso::TwoPartDate &t) noexcept {
    return std::upper_bound(tvec.begin(), tvec.end(), t);
  }

  /* @brief Reserve size/capacity for all instance's arrays 
   * @param[in] size Number of elements (per component) to be allocated
   */
  void reserve(int size) noexcept;

  /* @brief Add an EopRecord record to the instance; the function will add the
   *        new element in the correct position(s), ordered chornologically.
   * @param[in] rec The new EopRecord to be added to the instance.
   */
  void push_back(const EopRecord& rec) noexcept;

  /* @brief Constructor, given an optional size hint */
  EopLookUpTable(int size_hint=10) noexcept { reserve(size_hint); }

  /* @brief Destructor */
  ~EopLookUpTable() noexcept;

  /* @brief Return current number of elements/records in the instance */
  int size() const noexcept {return tvec.size();}

  /* @brief Overload operator()
   * @return The EopRecord at the index i
   */
  EopRecord operator()(int i) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < size());
#endif
    return EopRecord(
        {tvec[i], xp[i], yp[i], dut1[i], lod[i], dX[i], dY[i], xrt[i], yrt[i]});
  }

  /* @brief Get a reference to the epochs vector */
  const std::vector<dso::TwoPartDate> &epoch_vector() const noexcept {
    return tvec;
  }

  /* @brief Transform the instance's tvec vector from UTC to TT */
  void utc2tt() noexcept;

  /* @brief Compute (UT1-TAI)_n = (UT1-UTC)_n - ΔAT_n for n in tvec vector
   * 
   *  Transform the tabulated UT1-UTC (=dut1) values of the instance to 
   *  respective UT1-TAI values. The computed UT1-TAI are stored in the 
   *  returned vector (the instance is kept const), exactly corresponding to 
   *  the instance's tvec vector. Note that if you use this vector to 
   *  interpolate UT1-UTC values, you will need to add the ΔAT value for the 
   *  epoch of interpolation.
   *
   * @warning It is assumed that tvec contains UTC dates
   * 
   *  See Bradly et al, 2015, "Influence of ITRS/GCRS implementation for 
   *      astrodynamics: Coordinate transformations"
   */
  std::vector<double> dutToUt1Tai() const noexcept;
  
  /* @brief Compute (UT1-TAI)_n = (UT1-UTC)_n - ΔAT_n for n in tvec vector
   * 
   *  Transform the tabulated UT1-UTC (=dut1) values of the instance to 
   *  respective UT1-TAI values. The computed UT1-TAI are stored in the 
   *  the instance's tvec vector. Note that if you use this vector to 
   *  interpolate UT1-UTC values, you will need to add the ΔAT value for the 
   *  epoch of interpolation.
   *
   * @warning It is assumed that tvec contains UTC dates
   */
  void dutToUt1TaiInPlace() noexcept {dut1 = dutToUt1Tai();}

  /* @brief Evaluate the effects of zonal Earth tides on Earth rotation (i.e.
   *        ΔUT1 and ΔLOD) based on RG_ZONT2
   * We consider here zonal tidal variations with frequencies ranging from 
   * 5 days to 18.6 years.
   *
   * @param[in] dUt1 Effects of tides on dUT1 at the tabulated epochs [sec]
   * @param[in] dLod Effect on excess length of day (LOD) at the tabulated 
   *                 epochs [sec] (actually sec per day)
   * @param[in] tInUtc Signal that the t vector is in UTC time scale; else, it
   *            is assumed that t is in TT scale
   */
  void regularizedCorrections(std::vector<double> &dUt1,
                              std::vector<double> &dLod,
                              bool tInUtc) const noexcept;

  /* @brief Compute and remove the effect of zonal Earth tides on ΔUT1 and 
   *        ΔLOD, to get the so-called "regularized" corresponding values
   *        ΔUT1R and ΔLODR
   * 
   * According to Bradly et al, 2015: 
   * Prior to the interpolation of DUT1 and LOD, the tabulated values should 
   * be smoothed through regularization to enhance the interpolation accuracy. 
   * Regularization is the removal of zonal tidal variations with frequencies 
   * ranging from 5 days to 18.6 years. Regularized UT1 is typically denoted 
   * UT1R.
   *
   * Note that after regularization and interpolation, the zonal tide value 
   * should be added back at the time of interpolation.
   * 
   * All in all, this function will practically replace the arrays ΔUT1 and 
   * ΔLOD (i.e. dut and lod) with their regularized values ΔUT1-R and ΔLOD-R
   * 
   * @warning It is assumed that the calling instance's ΔUT1 and ΔLOD values 
   * are in units of [sec].
   * 
   * @param[in] tInUtc Signal that the t vector is in UTC time scale; else, it
   *            is assumed that t is in TT scale
   */
  void regularize(bool tInUtc) noexcept;

  /* @brief Perform Lagrangian interpolation of order 'order', to find the 
   *        EOP values at epoch 'mjd'. 
   * 
   * The function will fail if not enough data points are available for the 
   * interpolation (depending on the order parameter).
   *
   * @param[in] mjd   Date to interpolate at
   * @param[out] eopr Results of interpolation
   * @param[in] order Order of interpolation. The actual number of points used 
   *                  for the interpolation will be order+1. Use an odd number 
   *                  for order, e.g. order=5
   * @return Anything other than 0 denotes an error
   */
  int lagrange_interpolation(const dso::TwoPartDate &mjd, dso::EopRecord &eopr,
                             int order) const noexcept;


  /// @brief Interpolate and correct to get EOP/ERP parameters at given date
  ///
  /// See IERS Conventions 2010, Petit et al., Chapter 5.5
  ///
  /// @param[in] mjd Interpolation epoch given as mjd [TT] (assuming that
  ///            the calling instance's mjd array is in TT)
  /// @param[out] eopr EOP/ERP parameters at requested epoch, as computed by
  ///            Lagrangian interpolation and the removal of ocen-tide effects
  ///            and libration effects.
  /// @return Anything other than 0 is an error
  int interpolate(const dso::TwoPartDate &mjd, dso::EopRecord &eopr, int order,
                  bool regularized) const noexcept;
};/* EopLookUpTable */

namespace details {

/* An enumeration type to define EOP format/series */
enum class IersEopFormat {C0414, C0420};

/* @brief Read header and decide the format of the passed in C04 file 
 * @param[in] c04fn An IERS EOP file in either C04/14 or C04/20 format
 * @param[out] type Depending on the header, the type of the input file
 * @return Anything other than zero means that the file was not identified! In 
 *         this case, type should not be used.
 */
int choose_c04_series(const char *c04fn, IersEopFormat &type) noexcept;

/* @brief Given a C04/14 IERS eop file, parse it and store the EOP data for 
 *        the dates [start_mjd, end_mjd) to the passed-in eoptable instance.
 * 
 * Note that the EopLookUptable instance will be cleared before we start the 
 * parsing (all previous data stored by the instance, if any, will be lost).
 * The data are stored in the same units as in the input EOP file, that is
 * arcsec for xp, yp, dX and Dy and seconds for dut1 and lod. Dates are 
 * stored as UTC MJDs.
 *
 * Note that for C04/14 files, the EOP data for xrt and yrt are not available.
 * This function will set the corresponding elements to 0.
 *
 * See https://hpiers.obspm.fr/iers/eop/eopc04_14/eopc04_IAU2000.62-now and 
 * https://hpiers.obspm.fr/iers/eop/eopc04_14/updateC04.txt
 *
 * @param[in] c04fn An IERS C04/14 EOP file
 * @param[in] start_mjd Starting MJD to collect data for, in UTC (inclusive)
 * @param[in] end_mjd Last MJD of the interval to collect data for (exclusive, 
 *                    aka data will not be collected for this date) in UTC
 * @param[in] eoptable A EopLookUpTable where the parsed data will be stored 
 *            in. The instance will be cleared, hence any data stored in it 
 *            before the call, will be lost. Vectors of eoptable have the 
 *            following units:
 *            t         UTC as MJD at 0-hours in day
 *            xp, yp   (")
 *            dut      UT1 - UTC [sec]
 *            lod      Lenght of day offset LOD [sec]
 *            dx, dy   Celestial pole offsets dX, dY ["]
 *            xrt, yrt Pole coordinate rates ["/day], set to 0
 * @return At success 0 is returned. Else, something went wrong
 */
int parse_iers_C0414(const char *c04fn, dso::modified_julian_day start_mjd,
                        dso::modified_julian_day end_mjd,
                        dso::EopLookUpTable &eoptable) noexcept;

/* @brief Given a C04/20 IERS eop file, parse it and store the EOP data for 
 *        the dates [start_mjd, end_mjd) to the passed-in eoptable instance.
 * Note that the eoptable instance will be cleared before we start the 
 * parsing (all previous data stored by the instance, if any, will be lost).
 * The data are stored in the same units as in the input EOP file, that is
 * arcsec for xp, yp, dX and Dy and seconds for dut1 and lod. Dates are 
 * stored as UTC MJDs.
 * See also https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.txt and
 * https://hpiers.obspm.fr/iers/eop/eopc04_20/eopc04.1962-now
 * @param[in] c04fn An IERS C04/20 EOP file
 * @param[in] start_mjd Starting MJD to collect data for, in UTC (inclusive)
 * @param[in] end_mjd Last MJD of the interval to collect data for (exclusive, 
 *                    aka data will not be collected for this date) in UTC
 * @param[in] eoptable A EopLookUpTable where the parsed data will be stored 
 *            in. The instance will be cleared, hence any data stored in it 
 *            before the call, will be lost. Vectors of eoptable have the 
 *            following units:
 *            t         UTC as MJD at 0-hours in day
 *            xp, yp   (")
 *            dut      UT1 - UTC [sec]
 *            lod      Lenght of day offset LOD [sec]
 *            dx, dy   Celestial pole offsets dX, dY ["]
 *            xrt, yrt Pole coordinate rates ["/day]
 * @return At success 0 is returned. Else, something went wrong
 */
int parse_iers_C0420(const char *c04fn, dso::modified_julian_day start_mjd,
                        dso::modified_julian_day end_mjd,
                        dso::EopLookUpTable &eoptable) noexcept;
} /* namespace details */

inline int parse_iers_C04(const char *c04fn, dso::modified_julian_day start_mjd,
                          dso::modified_julian_day end_mjd,
                          dso::EopLookUpTable &eoptable) noexcept {
  details::IersEopFormat type;

  if (details::choose_c04_series(c04fn, type))
    return 1;

  if (type == details::IersEopFormat::C0420)
    return details::parse_iers_C0420(c04fn, start_mjd, end_mjd, eoptable);
  else
    return details::parse_iers_C0414(c04fn, start_mjd, end_mjd, eoptable);
}
} /* namespace dso */

#endif
