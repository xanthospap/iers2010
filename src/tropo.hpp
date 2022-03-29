#ifndef __TROPO_MODELS_ALTERNATIVES_HPP__
#define __TROPO_MODELS_ALTERNATIVES_HPP__

///
/// @file tropo.hpp
/// Declerations for a number of functions that are not included in the
/// IERS2010 standards but are used in satellite techniques to handle
/// path delay caused by troposphere.
/// All of the declerations are nested inside the 'dso' namespace (and not 
/// 'iers2010').
///

#include "datetime/dtcalendar.hpp"
#ifdef DEBUG
#include <cstdio>
#endif

namespace dso {

/// @brief Zenith Wet Delay by Aske & Nordius, 1987
/// This function determines the zenith wet delay based on the
/// equation 22 by Aske and Nordius (1987)
/// Reference:
/// Askne and Nordius, Estimation of tropospheric delay for microwaves from
/// surface weather data, Radio Science, Vol 22(3): 379-386, 1987.
/// Translated from the MATLAB source code at TU Vienna:
/// https://vmf.geo.tuwien.ac.at/codes/asknewet.m
///
/// @param[in] e   water vapor pressure in hPa
/// @param[in] Tm  mean temperature in Kelvin
/// @param[in] lambda water vapor lapse rate (see definition in Askne and
///            Nordius 1987)
/// @return zwd:  zenith wet delay in meters
double asknewet(double e, double Tm, double lambda) noexcept;

/// @brief Compute zenith hydrostatic delay using refined Saastamoinen
/// This subroutine determines the zenith hydrostatic delay based on the
/// equation by Saastamoinen (1972) as refined by Davis et al. (1985).
/// Translated from MATLAB source, originaly found at:
/// https://vmf.geo.tuwien.ac.at/codes/saasthyd.m
/// Saastamoinen, J., Atmospheric correction for the troposphere and
/// stratosphere in radio ranging of satellites. The use of artificial
/// satellites for geodesy, Geophys. Monogr. Ser. 15, Amer. Geophys. Union,
/// pp. 274-251, 1972.
/// Davis, J.L, T.A. Herring, I.I. Shapiro, A.E.E. Rogers, and G. Elgered,
/// Geodesy by Radio Interferometry: Effects of Atmospheric Modeling Errors
/// on Estimates of Baseline Length, Radio Science, Vol. 20, No. 6,
/// pp. 1593-1607, 1985.
/// param[in] p     pressure in hPa
/// param[in] dlat  ellipsoidal latitude in radians
/// param[in] dlon  longitude in radians
/// param[in] hell  ellipsoidal height in m
/// return zhd: zenith hydrostatic delay in meter
double saasthyd(double p, double dlat, double hell) noexcept;
int vmf3(double ah, double aw, dso::datetime<dso::nanoseconds> &t, double lat,
         double lon, double zd, double &mfh, double &mfw) noexcept;

/// @brief gpt3 namespace holds declerations of gpt3 implementation
namespace gpt3 {

/// @enum We can have two different grid files for gpt3, based on the step
///       size of the grid; one is 1x1 deg. and one is 5x5 deg.
enum class Gpt3Grid : char { grid1x1, grid5x5 };

/// @class gpt3_grid_attributes Holds attributes for each individual grid type
template <Gpt3Grid G> struct gpt3_grid_attributes {};

/// @class gpt3_grid_attributes Holds attributes for the 1x1 grid type
template <> struct gpt3_grid_attributes<Gpt3Grid::grid1x1> {
  template <typename T> static constexpr T grid_tick() { return T(1e0); }
  static constexpr unsigned num_lines = 64801 - 1;
};

/// @class gpt3_grid_attributes Holds attributes for the 5x5 grid type
template <> struct gpt3_grid_attributes<Gpt3Grid::grid5x5> {
  template <typename T> static constexpr T grid_tick() { return T(5e0); }
  static constexpr unsigned num_lines = 2593 - 1;
};

/// @class gpt3_grid
/// A class to hold the data arrays parsed from the GPT3 grid file(s) and
/// needed to perform gpt3-related computations.
/// @warning At least for the 1x1 grid, this class is too large to be held at
///          stack (causes stack overflow on some occasions), hence it is 
///          allocated on the heap.
/// @note At least for now, this class is not copyable or movable; just pass
///       it around via pointers (or references).
struct gpt3_grid {

  /// @brief Initialize and allocate given the number of rows of the 
  ///        corresponding grid file.
  /// @param[in] rows Number of rows of the corresponding grid file
  gpt3_grid(int rows = 0) {
    if (rows)
      allocate(rows);
  };

  /// @brief Destructor; frees allocated memory
  ~gpt3_grid() noexcept {
    if (size)
      dealloc();
  }

  /// @brief Copy not allowed
  gpt3_grid(const gpt3_grid &) noexcept = delete;
  
  /// @brief Move not allowed
  gpt3_grid(gpt3_grid &&) noexcept = delete;
  
  /// @brief Copy-assignment not allowed
  gpt3_grid &operator=(const gpt3_grid &) noexcept = delete;
  
  /// @brief Move-assignment not allowed
  gpt3_grid &operator=(gpt3_grid &&) noexcept = delete;

  /// @brief Perform allocations
  /// @param[in] num_rows Number of rows for each instance array (aka the
  ///            number of rows of the corresponding grid file)
  void allocate(unsigned num_rows);

  /// @brief Free allocated memory
  void dealloc() noexcept;

  unsigned size = 0;
  double **p_grid = nullptr;
  double **T_grid = nullptr;
  double **Q_grid = nullptr;
  double **dT_grid = nullptr;
  double *u_grid = nullptr;
  double *Hs_grid = nullptr;
  double **ah_grid = nullptr;
  double **aw_grid = nullptr;
  double **la_grid = nullptr;
  double **Tm_grid = nullptr;
  double **Gn_h_grid = nullptr;
  double **Ge_h_grid = nullptr;
  double **Gn_w_grid = nullptr;
  double **Ge_w_grid = nullptr;
};

/// @brief Parse a GPT3 grid file (either 1x1 or 5x5 degrees)
/// This functionwill parse a GPT3 grid file, either 1x1 or 5x5 degrees, and
/// store the values in the corresponding arrays of the passed-in grid
/// instance.
/// The grid files can be found at: https://vmf.geo.tuwien.ac.at/codes/
/// and are respectively:
/// - gpt3_5.grd and
/// - gpt3_1.grd
/// The format of these files is specific and the function expects that the
/// files adhere to it. It also expexts to find a header line (at the top of 
/// the input file).
/// @param[in] gridfn The filename of the GPT3 grid file to be read
/// @param[out] grid  A gpt3_grid instance to hold parsed values. The gpt3_grid
///            instance should have the correct size depending on the file to
///            be read, aka the instance should have already been constructed
///            with the correct size before passed in.
/// @return Anything other than 0 denotes an error
int parse_gpt3_grid(const char *gridfn, gpt3_grid *grid) noexcept;
} // gpt3

/// @class gpt3_result A structure to hold GPT3 details returned by the
///        corresponding function.
struct gpt3_result {
  double p;    ///< pressure in hPa
  double T;    ///< temperature in degrees Celsius
  double dT;   ///< temperature lapse rate in degrees per km
  double Tm;   ///< mean temperature weighted with the water vapor in degrees
               ///< Kelvin
  double e;    ///< water vapour pressure in hPa
  double ah;   ///< hydrostatic mapping function coefficient(VMF3)
  double aw;   ///< wet mapping function coefficient(VMF3)
  double la;   ///< water vapour decrease factor
  double undu; ///< geoid undulation in m
  double Gn_h; ///< hydrostatic north gradient in m
  double Ge_h; ///< hydrostatic east gradient in m
  double Gn_w; ///< wet north gradient in m
  double Ge_w; ///< wet east gradient in m
};// gpt3_result

/// @class vmf3_hw A structure to hold VMF3 details returned by the
///        corresponding function.
struct vmf3_hw {
  double mfh; ///< Hydrostatic mapping function
  double mfw; ///< Wet mapping function
}; // vmf3_result

int gpt3_fast(double fractional_doy, const double *lat,
                   const double *lon, const double *hell, int num_stations,
                   int it, const gpt3::gpt3_grid *gridNxN,
                   dso::gpt3_result *g3out) noexcept;
/// @brief Implement GPT3 for a number of stations
/// This subroutine determines pressure, temperature, temperature lapse rate, 
/// mean temperature of the water vapor, water vapour pressure, hydrostatic 
/// and wet mapping function coefficients ah and aw, water vapour decrease
/// factor, geoid undulation and empirical tropospheric gradients for 
/// specific sites near the earth's surface. All output values are valid for
/// the specified ellipsoidal height hell.
/// GPT3_5 is based on a 1x1 or 5x5 external grid file ('gpt3_[51].grd') with
/// mean values as well as sine and cosine amplitudes for the annual and
/// semiannual variation of the coefficients.
/// Note that this version is 'agnostic' as to what grid file is used; it takes
/// as an input parameter a gpt3_grid instance, which should already hold all 
/// the needed values. The function only know how many rows this instance has 
/// (and hence can conclude if it is the 1x1 or the 5x5 version).
/// Translated from MATLAB source, originaly found at:
/// https://vmf.geo.tuwien.ac.at/codes/gpt3_[51]_fast.m
/// Reference:
/// D. Landskron, J. Bohm (2018), VMF3/GPT3: Refined Discrete and Empirical 
/// Troposphere Mapping Functions, J Geod (2018) 92: 349., 
/// doi: 10.1007/s00190-017-1066-2. 
/// Download at: 
/// https://link.springer.com/content/pdf/10.1007%2Fs00190-017-1066-2.pdf
/// Translated from the gpt3_[15]_fast.m MATLAB source code, found at:
/// https://vmf.geo.tuwien.ac.at/codes/ by TU Vienna.
///
/// @param[in] t     the date to perform computations for
/// @param[in] lat   ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
/// @param[in] lon   longitude in radians [-pi:pi] or [0:2pi] (vector)
/// @param[in] h_ell ellipsoidal height in m (vector)
/// @param[in] num_stations Number of stations, aka the size of the lat, lon
///                  and hell arrays
/// @param[in] it    1: no time variation but static quantities
///                  0: with time variation (annual and semiannual terms)
/// @param[in] gridNxN A gpt3_grid instance, holding the needed values to 
///                  perform the computation. This instance should already
///                  hold the needed values, aka the corresponding grid fil
///                  should have been parsed and values stored in gridNxN.
/// @param[out] g3out An instance of gpt3_result where all computed values are
///                  stored at.
/// @return Anything other than 0 denotes an error
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
int gpt3_fast(const dso::datetime<S> &t, const double *lat,
              const double *lon, const double *hell, int num_stations, int it,
              const gpt3::gpt3_grid *gridNxN, gpt3_result *g3out) noexcept {
 const auto yrdoy = t.mjd().to_ydoy();
 double fraction = t.sec().fractional_days();
 fraction += static_cast<double>(yrdoy.__doy.as_underlying_type());
 return gpt3_fast(fraction, lat, lon, hell, num_stations, it, gridNxN, g3out);
}

/// @brief Implement GPT3 for a number of stations
/// This subroutine determines pressure, temperature, temperature lapse rate, 
/// mean temperature of the water vapor, water vapour pressure, hydrostatic 
/// and wet mapping function coefficients ah and aw, water vapour decrease
/// factor, geoid undulation and empirical tropospheric gradients for 
/// specific sites near the earth's surface. All output values are valid for
/// the specified ellipsoidal height hell.
/// GPT3_5 is based on a 1x1 or 5x5 external grid file ('gpt3_[51].grd') with
/// mean values as well as sine and cosine amplitudes for the annual and
/// semiannual variation of the coefficients.
/// Translated from MATLAB source, originaly found at:
/// https://vmf.geo.tuwien.ac.at/codes/gpt3_[51]_fast.m
/// Reference:
/// D. Landskron, J. Bohm (2018), VMF3/GPT3: Refined Discrete and Empirical 
/// Troposphere Mapping Functions, J Geod (2018) 92: 349., 
/// doi: 10.1007/s00190-017-1066-2. 
/// Download at: 
/// https://link.springer.com/content/pdf/10.1007%2Fs00190-017-1066-2.pdf
/// Translated from the gpt3_[15]_fast.m MATLAB source code, found at:
/// https://vmf.geo.tuwien.ac.at/codes/ by TU Vienna.
///
/// @param[in] t     the date to perform computations for
/// @param[in] lat   ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
/// @param[in] lon   longitude in radians [-pi:pi] or [0:2pi] (vector)
/// @param[in] h_ell ellipsoidal height in m (vector)
/// @param[in] num_stations Number of stations, aka the size of the lat, lon
///                  and hell arrays
/// @param[in] it    1: no time variation but static quantities
///                  0: with time variation (annual and semiannual terms)
/// @param[in]  grid_file The filename of the grid file to be used, either
///                  'gpt3_5.grd' or 'gpt3_1.grd' depending on the resolution
///                  wanted
/// @param[out] g3out An instance of gpt3_result where all computed values are
///                  stored at.
/// @return Anything other than 0 denotes an error
int gpt3_fast(const dso::datetime<dso::nanoseconds> &t, const double *lat,
              const double *lon, const double *hell, int num_stations, int it,
              const char *grid_file, gpt3_result *g3out,
              int &grid_step) noexcept;

/// @brief Determine the VMF3 hydrostatic and wet mapping factors
/// This function determines the VMF3 hydrostatic and wet mapping factors.
/// The a coefficients have to be inserted from discrete data, while the b
/// and c coefficients are of empirical nature containing a geographical 
/// and temporal dependence, represented in spherical harmonics. The 
/// spherical harmonics coefficients are developed to degree and order 12 and 
/// are based on a 5x5 grid containing ray-tracing data from 2001-2010.
/// This function needs the ah and aw values (aka hydrostatic mf coefficient a
/// and wet mf coefficient a) which can be obtained from:
/// http://vmf.geo.tuwien.ac.at/trop_products/ but in this implementation they
/// are read off from an dso::gpt3_result instance passed in.
/// Translated from the vmf3.m MATLAB source code, found at:
/// https://vmf.geo.tuwien.ac.at/codes/ by TU Vienna.
///
/// @param[in] gptres An array of dso::gpt3_result where values for 
///            ah (hydrostatic mf coefficient a) and aw (wet mf coefficient a)
///            values are stored
/// @param[in] The epoch for which we want the computation
/// @param[in] lat latitude (radians)
/// @param[in] lon longitude (radians)
/// @param[in] zd zenith distance (radians)
/// @param[in] num_sta number of stations involved, aka the size of the
///            gptres, vlat, vlon, vmfres and vzd arrays
/// @param[out] vmfres results array; each element includes the values for
///            - mfh: hydrostatic mapping factor and
///            - mfw: wet mapping factor
int vmf3(const dso::gpt3_result *gptres, dso::datetime<dso::nanoseconds> &t,
         const double *vlat, const double *vlon, const double *vzd,
         dso::vmf3_hw *vmfres, int num_sta) noexcept;

} // dso
#endif
