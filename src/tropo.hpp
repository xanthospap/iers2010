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
#include <stdexcept>
#include <vector>
#include <array>

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
/// @param[in] e   water vapor pressure [hPa]
/// @param[in] Tm  mean temperature [Kelvin]
/// @param[in] lambda water vapor lapse rate (see definition in Askne and
///            Nordius 1987)
/// @return zwd: zenith wet delay [meters]
double asknewet(double e, double Tm, double lambda) noexcept;

/// @brief Compute zenith hydrostatic delay using refined Saastamoinen model.
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
/// @param[in] p     pressure in [hPa]
/// @param[in] dlat  ellipsoidal latitude [radians]
/// @param[in] dlon  longitude [radians]
/// @param[in] hell  ellipsoidal height [meters]
/// @return zhd: zenith hydrostatic delay [meters]
double saasthyd(double p, double dlat, double hell) noexcept;

int vmf3(double ah, double aw, const dso::datetime<dso::nanoseconds> &t,
         double lat, double lon, double zd, double &mfh, double &mfw) noexcept;
int vmf3_ht(double ah, double aw, const dso::datetime<dso::nanoseconds> &t,
            double lat, double lon, double h_ell, double zd, double &mfh,
            double &mfw) noexcept;

/// @brief gpt3 namespace holds declerations of gpt3 implementation
namespace gpt3 {

/// @enum We can have two different grid files for gpt3, based on the step
///       size of the grid; one is 1x1 deg. and one is 5x5 deg.
enum class Gpt3GridResolution : char { grid1x1, grid5x5 };

/// @class gpt3_grid_attributes Holds attributes for each individual grid type
template <Gpt3GridResolution G> struct gpt3_grid_attributes {};

/// @class gpt3_grid_attributes Holds attributes for the 1x1 grid type
template <> struct gpt3_grid_attributes<Gpt3GridResolution::grid1x1> {
  template <typename T> static constexpr T grid_tick() { return T(1e0); }
  static constexpr unsigned num_lines = 64801 - 1;
};

/// @class gpt3_grid_attributes Holds attributes for the 5x5 grid type
template <> struct gpt3_grid_attributes<Gpt3GridResolution::grid5x5> {
  template <typename T> static constexpr T grid_tick() { return T(5e0); }
  static constexpr unsigned num_lines = 2593 - 1;
};
}// gpt3

/// @class Gpt3Grid
/// A class to hold the data arrays parsed from the GPT3 grid file(s) and
/// needed to perform gpt3-related computations.
/// @warning At least for the 1x1 grid, this class is too large to be held at
///          stack (causes stack overflow on some occasions), hence it is 
///          allocated on the heap.
/// @note At least for now, this class is not copyable or movable; just pass
///       it around via pointers (or references).
struct Gpt3Grid {
  double *memPool{nullptr};
  int num_rows{0};
  unsigned int offset{0};
  static constexpr const int cols = 5;
  // arrays **MUST** be in the following order:
  // p_offset
  // t_offset
  // q_offset
  // dt_offset
  // ah_offset
  // aw_offset
  // la_offset
  // tm_offset
  // gn_h_offset
  // ge_h_offset
  // gn_w_offset
  // ge_w_offset
  // u_offset
  // hs_offset

/// @brief Parse a GPT3 grid file (either 1x1 or 5x5 degrees)
/// This function will parse a GPT3 grid file, either 1x1 or 5x5 degrees, and
/// store the values in the corresponding arrays of the calling grid instance.
/// The grid files can be found at: https://vmf.geo.tuwien.ac.at/codes/
/// and are respectively:
/// - gpt3_5.grd and
/// - gpt3_1.grd
/// The format of these files is specific and the function expects that the
/// files adhere to it. It also expexts to find a header line (at the top of 
/// the input file).
/// Note that the function will change the calling instance in two ways:
/// 1. set's it properties (member variables and memmory) tyo the right values
/// 2. parse the grid file and fill in the instance's matrices/arrays
///
/// @param[in] gridfn The filename of the GPT3 grid file to be read
/// @return Anything other than 0 denotes an error
  int parse_grid(const char *fn) noexcept;

  Gpt3Grid() noexcept {};
  Gpt3Grid(const char *fn) {
    if (parse_grid(fn))
      throw std::runtime_error("Failed to create Gpt3Grid instance");
  }
  
  int allocate(double resolution);

  void deallocate() noexcept {
    if (memPool) delete[] memPool;
    memPool = nullptr;
    num_rows = 0;
    offset = 0;
  }

  ~Gpt3Grid() noexcept {
    if (memPool) delete[] memPool;
  }

  double *p_grid_row(int row) noexcept {
    return memPool + 0 * offset + cols * row;
  }
  const double *p_grid_row(int row) const noexcept {
    return memPool + 0 * offset + cols * row;
  }
  double p_grid(int row, int col) const noexcept {
    return p_grid_row(row)[col];
  }

  double *t_grid_row(int row) noexcept {
    return memPool + 1 * offset + cols * row;
  }
  const double *t_grid_row(int row) const noexcept {
    return memPool + 1 * offset + cols * row;
  }
  double t_grid(int row, int col) const noexcept {
    return t_grid_row(row)[col];
  }

  double *q_grid_row(int row) noexcept {
    return memPool + 2 * offset + cols * row;
  }
  const double *q_grid_row(int row) const noexcept {
    return memPool + 2 * offset + cols * row;
  }
  double q_grid(int row, int col) const noexcept {
    return q_grid_row(row)[col];
  }

  double *dt_grid_row(int row) noexcept {
    return memPool + 3 * offset + cols * row;
  }
  const double *dt_grid_row(int row) const noexcept {
    return memPool + 3 * offset + cols * row;
  }
  double dt_grid(int row, int col) const noexcept {
    return dt_grid_row(row)[col];
  }

  double *ah_grid_row(int row) noexcept {
    return memPool + 4 * offset + cols * row;
  }
  const double *ah_grid_row(int row) const noexcept {
    return memPool + 4 * offset + cols * row;
  }
  double ah_grid(int row, int col) const noexcept {
    return ah_grid_row(row)[col];
  }

  double *aw_grid_row(int row) noexcept {
    return memPool + 5 * offset + cols * row;
  }
  const double *aw_grid_row(int row) const noexcept {
    return memPool + 5 * offset + cols * row;
  }
  double aw_grid(int row, int col) const noexcept {
    return aw_grid_row(row)[col];
  }

  double *la_grid_row(int row) noexcept {
    return memPool + 6 * offset + cols * row;
  }
  const double *la_grid_row(int row) const noexcept {
    return memPool + 6 * offset + cols * row;
  }
  double la_grid(int row, int col) const noexcept {
    return la_grid_row(row)[col];
  }

  double *tm_grid_row(int row) noexcept {
    return memPool + 7 * offset + cols * row;
  }
  const double *tm_grid_row(int row) const noexcept {
    return memPool + 7 * offset + cols * row;
  }
  double tm_grid(int row, int col) const noexcept {
    return tm_grid_row(row)[col];
  }

  double *gn_h_grid_row(int row) noexcept {
    return memPool + 8 * offset + cols * row;
  }
  const double *gn_h_grid_row(int row) const noexcept {
    return memPool + 8 * offset + cols * row;
  }
  double gn_h_grid(int row, int col) const noexcept {
    return gn_h_grid_row(row)[col];
  }

  double *ge_h_grid_row(int row) noexcept {
    return memPool + 9 * offset + cols * row;
  }
  const double *ge_h_grid_row(int row) const noexcept {
    return memPool + 9 * offset + cols * row;
  }
  double ge_h_grid(int row, int col) const noexcept {
    return ge_h_grid_row(row)[col];
  }

  double *gn_w_grid_row(int row) noexcept {
    return memPool + 10 * offset + cols * row;
  }
  const double *gn_w_grid_row(int row) const noexcept {
    return memPool + 10 * offset + cols * row;
  }
  double gn_w_grid(int row, int col) const noexcept {
    return gn_w_grid_row(row)[col];
  }

  double *ge_w_grid_row(int row) noexcept {
    return memPool + 11 * offset + cols * row;
  }
  const double *ge_w_grid_row(int row) const noexcept {
    return memPool + 11 * offset + cols * row;
  }
  double ge_w_grid(int row, int col) const noexcept {
    return ge_w_grid_row(row)[col];
  }

  double *u_grid_row() noexcept { return memPool + 12 * offset; }
  const double *u_grid_row() const noexcept {
    return memPool + 12 * offset;
  }
  double u_grid(int col) const noexcept {
    return u_grid_row()[col];
  }

  double *hs_grid_row() noexcept {
    return memPool + 12 * offset + num_rows;
  }
  const double *hs_grid_row() const noexcept {
    return memPool + 12 * offset + num_rows;
  }
  double hs_grid(int col) const noexcept {
    return hs_grid_row()[col];
  }

  /// @brief Copy not allowed
  Gpt3Grid(const Gpt3Grid &) noexcept = delete;
  
  /// @brief Move not allowed
  Gpt3Grid(Gpt3Grid &&) noexcept = delete;
  
  /// @brief Copy-assignment not allowed
  Gpt3Grid &operator=(const Gpt3Grid &) noexcept = delete;
  
  /// @brief Move-assignment not allowed
  Gpt3Grid &operator=(Gpt3Grid &&) noexcept = delete;

}; //Gpt3Grid

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
/// @param[in] fractional_doy the date as fractional day of year to perform 
///            computations for
/// @param[in] ellipsoidal An array of 3-dimensional arrays containing 
///            ellipsoidal coordinates for num_stations stations. These arrays
///            should hold longitude/latitude/height as:
///            * ellipsoidal latitude in range (-pi/2:+pi/2), [radians]
///            * longitude in range (-pi:pi) or (0:2pi), [radians]
///            * ellipsoidal height [meters]
/// @param[in] num_stations Number of stations, aka the size of ellipsoidal
/// @param[in] it    1: no time variation but static quantities
///                  0: with time variation (annual and semiannual terms)
/// @param[in] grid A Gpt3Grid instance, holding the needed values to 
///                  perform the computation. This instance should already
///                  hold the needed values, aka the corresponding grid file
///                  should have been parsed and values stored in gridNxN.
/// @param[out] g3out An instance of gpt3_result where all computed values are
///                  stored at. Should be at least of size num_stations.
/// @return Anything other than 0 denotes an error
int gpt3_fast(double fractional_doy,
              const std::vector<std::array<double, 3>> &ellipsoidal,
              int it, const Gpt3Grid &grid,
              std::vector<dso::gpt3_result> &g3out) noexcept;

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
/// @param[in] ellipsoidal An array of 3-dimensional arrays containing 
///            ellipsoidal coordinates for num_stations stations. These arrays
///            should hold longitude/latitude/height as:
///            * ellipsoidal latitude in range (-pi/2:+pi/2), [radians]
///            * longitude in range (-pi:pi) or (0:2pi), [radians]
///            * ellipsoidal height [meters]
/// @param[in] num_stations Number of stations, aka the size of ellipsoidal
/// @param[in] it    1: no time variation but static quantities
///                  0: with time variation (annual and semiannual terms)
/// @param[in] gridNxN A gpt3_grid instance, holding the needed values to 
///                  perform the computation. This instance should already
///                  hold the needed values, aka the corresponding grid file
///                  should have been parsed and values stored in gridNxN.
/// @param[out] g3out An instance of gpt3_result where all computed values are
///                  stored at. Should be at least of size num_stations.
/// @return Anything other than 0 denotes an error
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt S>
#else
template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
int gpt3_fast(const dso::datetime<S> &t, const std::vector<std::array<double,3>> &ellipsoidal,
              int it, const Gpt3Grid &grid,
              std::vector<dso::gpt3_result> &g3out) noexcept {
 
 const auto yrdoy = t.mjd().to_ydoy();
 double fraction = t.sec().fractional_days();
 fraction += static_cast<double>(yrdoy.__doy.as_underlying_type());
 
 return gpt3_fast(fraction, ellipsoidal, it, grid, g3out);
}


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
/*int vmf3(const dso::gpt3_result *gptres, dso::datetime<dso::nanoseconds> &t,
         const double *vlat, const double *vlon, const double *vzd,
         dso::vmf3_hw *vmfres, int num_sta) noexcept;*/

} // dso
#endif
