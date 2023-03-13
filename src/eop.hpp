#ifndef __DSO_IERS_PRODUCTS_EOP_HPP__
#define __DSO_IERS_PRODUCTS_EOP_HPP__

#include "datetime/dtcalendar.hpp"
#include "iers2010.hpp"
#include <vector>

namespace dso {

/// @bref A structure to hold EOP records for a single epoch
struct EopRecord {
  dso::TwoPartDate mjd; // MJD
  double xp, yp;  ///< Pole coordinates x, y (")
  double dut;     ///< UT1 - UTC [sec]
  double lod;     ///< Lenght of day offset LOD [sec]
  double dx, dy;  ///< Celestial pole offsets dX, dY ["]
  double xrt, yrt; ///< Pole coordinate rates ["/day], only available in C04/20
  
  /// @brief Angular velocity of Earth in [rad/sec], including LOD variation
  /// @return Angular velocity of Earth (Ω) in [rad/sec]
  double omega() const noexcept {
    // see https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html
    const double OmPicoRadSec =
        iers2010::OmegaEarth * 1e12 - 0.843994809e0 * (lod * 1e3);
    return OmPicoRadSec * 1e-12;
  }
};// EopRecord

/// @brief A class to hold EOP information ordered chronologically.
class EopLookUpTable {
  std::vector<dso::TwoPartDate> t; // UTC or TT (MJD)
  std::vector<double> xp, yp, dut1, dX, dY, lod, xrt, yrt;

public:
  // @brief Clear all instance's arrays
  void clear() noexcept;

  void printxp() const noexcept { 
    printf("Xp Series:\n");
    for (const auto x : xp) printf("%+9.6f, ", x);
    printf("\n");
  }
  
  // @brief Reserve size/capacity for all instance's arrays
  void reserve(int size) noexcept;

  void push_back(const EopRecord& rec) noexcept;

  EopLookUpTable(int size_hint=10) noexcept { reserve(size_hint); }
  ~EopLookUpTable() noexcept;
  int size() const noexcept {return t.size();}

  /// @brief Transform the instance's t vector from UTC to TT
  void utc2tt() noexcept;

  /// @brief Compute (UT1-TAI)_n = (UT1-UTC)_n - ΔAT_n for n in t vector
  ///
  /// Transform the tabulated UT1-UTC (=dut1) values of the instance to 
  /// respective UT1-TAI values. The computed UT1-TAI are stored in the 
  /// returned vector (the instance is kept const), exactly corresponding to 
  /// the instance's t vector. Note that if you use this vector to interpolate
  /// UT1-UTC values, you will need to add the ΔAT value for the epoch of 
  /// interpolation.
  ///
  /// See Bradly et al, 2015, "Influence of ITRS/GCRS implementation for 
  ///     astrodynamics: Coordinate transformations"
  std::vector<double> dutToUt1Tai() const noexcept;
  void dutToUt1TaiInPlace() noexcept {dut1 = dutToUt1Tai();}

  /// @brief Evaluate the effects of zonal Earth tides on Earth rotation (i.e.
  ///        ΔUT1 and ΔLOD) based on RG_ZONT2
  /// We consider here zonal tidal variations with frequencies ranging from 
  /// 5 days to 18.6 years.
  /// @param[in] dUt1 Effects of tides on dUT1 at the tabulated epochs [sec]
  /// @param[in] dLod Effect on excess length of day (LOD) at the tabulated 
  ///                 epochs [sec] (actually sec per day)
  /// @param[in] tInUtc Signal that the t vector is in UTC time scale; else, it
  ///            is assumed that t is in TT scale
  void regularizedCorrections(std::vector<double> &dUt1,
                              std::vector<double> &dLod,
                              bool tInUtc) const noexcept;

  /// @brief Compute and remove the effect of zonal Earth tides on ΔUT1 and 
  ///        ΔLOD, to get the so-called "regularized" corresponding values
  ///        ΔUT1R and ΔLODR
  /// According to Bradly et al, 2015: 
  /// Prior to the interpolation of DUT1 and LOD, the tabulated values should 
  /// be smoothed through regularization to enhance the interpolation accuracy. 
  /// Regularization is the removal of zonal tidal variations with frequencies 
  /// ranging from 5 days to 18.6 years. Regularized UT1 is typically denoted 
  /// UT1R.
  /// Note that after regularization and interpolation, the zonal tide value 
  /// should be added back at the time of interpolation.
  /// All in all, this function will practically replace the arrays 
  /// ΔUT1 and ΔLOD (i.e. dut and lod) with their regularized values ΔUT1-R 
  /// and ΔLOD-R
  /// @warning It is assumed that the calling instance's ΔUT1 and ΔLOD values 
  /// are in units of [sec].
  /// @param[in] tInUtc Signal that the t vector is in UTC time scale; else, it
  ///            is assumed that t is in TT scale
  void regularize(bool tInUtc) noexcept;

  /// @brief Perform Lagrangian interpolation of order 'order',  to find the 
  ///        EOP values at epoch 'mjd'
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
};// EopLookUpTable

/// @brief Given a C04/14 IERS eop file, parse it and store the EOP data for 
///        the dates [start_mjd, end_mjd) to the passed-in eoptable instance.
/// Note that the eoptable instance will be cleared before we start the 
/// parsing (all previous data stored by the instance, if any, will be lost).
/// The data are stored in the same units as in the input EOP file, that is
/// arcsec for xp, yp, dX and Dy and seconds for dut1 and lod. Dates are 
/// stored as UTC MJDs.
/// Note that for C04/14 files, the EOP data for xrt and yrt are not available.
/// This function will set the corresponding elements to 0.
/// See https://hpiers.obspm.fr/iers/eop/eopc04_14/eopc04_IAU2000.62-now and 
/// https://hpiers.obspm.fr/iers/eop/eopc04_14/updateC04.txt
/// @param[in] c04fn An IERS C04/14 EOP file
/// @param[in] start_mjd Starting MJD to collect data for, in UTC (inclusive)
/// @param[in] end_mjd Last MJD of the interval to collect data for (exclusive, 
///                    aka data will not be collected for this date) in UTC
/// @param[in] eoptable A EopLookUpTable where the parsed data will be stored 
///            in. The instance will be cleared, hence any data stored in it 
///            before the call, will be lost. Vectors of eoptable have the 
///            following units:
///            t         UTC as MJD at 0-hours in day
///            xp, yp   (")
///            dut      UT1 - UTC [sec]
///            lod      Lenght of day offset LOD [sec]
///            dx, dy   Celestial pole offsets dX, dY ["]
///            xrt, yrt Pole coordinate rates ["/day], set to 0
/// @return At success 0 is returned. Else, something went wrong
int parse_iers_C0414(const char *c04fn, dso::modified_julian_day start_mjd,
                        dso::modified_julian_day end_mjd,
                        dso::EopLookUpTable &eoptable) noexcept;

/// @brief Given a C04/20 IERS eop file, parse it and store the EOP data for 
///        the dates [start_mjd, end_mjd) to the passed-in eoptable instance.
/// Note that the eoptable instance will be cleared before we start the 
/// parsing (all previous data stored by the instance, if any, will be lost).
/// The data are stored in the same units as in the input EOP file, that is
/// arcsec for xp, yp, dX and Dy and seconds for dut1 and lod. Dates are 
/// stored as UTC MJDs.
/// See also https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.txt and
/// https://hpiers.obspm.fr/iers/eop/eopc04_20/eopc04.1962-now
/// @param[in] c04fn An IERS C04/20 EOP file
/// @param[in] start_mjd Starting MJD to collect data for, in UTC (inclusive)
/// @param[in] end_mjd Last MJD of the interval to collect data for (exclusive, 
///                    aka data will not be collected for this date) in UTC
/// @param[in] eoptable A EopLookUpTable where the parsed data will be stored 
///            in. The instance will be cleared, hence any data stored in it 
///            before the call, will be lost. Vectors of eoptable have the 
///            following units:
///            t         UTC as MJD at 0-hours in day
///            xp, yp   (")
///            dut      UT1 - UTC [sec]
///            lod      Lenght of day offset LOD [sec]
///            dx, dy   Celestial pole offsets dX, dY ["]
///            xrt, yrt Pole coordinate rates ["/day]
/// @return At success 0 is returned. Else, something went wrong
int parse_iers_C0420(const char *c04fn, dso::modified_julian_day start_mjd,
                        dso::modified_julian_day end_mjd,
                        dso::EopLookUpTable &eoptable) noexcept;
}// namespace dso

#endif
