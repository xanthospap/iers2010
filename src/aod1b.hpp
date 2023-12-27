/** @file
 * Utilities to read/parse AOD1B files
 */

#ifndef __DSO_AOD1B_IPARSER_GPP__
#define __DSO_AOD1B_IPARSER_GPP__

#include "datetime/calendar.hpp"
#include <cstring>

namespace dso {

/** Different sets of Stokes coefficients included in AOD1B products */
enum class AOD1BCoefficientType { ATM, OCN, GLO, OBA };

class Aod1bIn {
  /* from offset=0 with size = 19 + '\0' */
  static constexpr int agency_offset = 0;
  /* from offset=20 with size = 19 + '\0' */
  static constexpr int satellite_offset = 20;
  /* from offset=40 with size = 19 + '\0' */
  static constexpr int sensor_offset = 40;
  /* from offset=60 with size = 3 + '\0' */
  static constexpr int process_level_offset = 60;
  /* from offset=64 with size = 3 + '\0' */
  static constexpr int pressure_type_offset = 64;
  /* max number of chars in charArena */
  static constexpr int MaxArenaChars = 80;

private:
  /* filename */
  std::string mfn;
  /* buffer to store all character/string information in an AOD1B header */
  char charArena[MaxArenaChars];
  /* file type -- should be 999 */
  int mfile_type;
  /* file format: 0=binary 1=ascii */
  int mfile_format;
  /* number of header records */
  int mnum_header_records;
  /* number of data records */
  int mnum_data_records;
  /* maximum degree for harmonic coeffs */
  int mmax_degree;
  /* coefficient errors (yes/no) */
  int mcoeff_errors;
  /* coeff. normalized (yes/no) */
  int mcoeff_normalized;
  /* constant gm [m^3/s^2] */
  double mGM;
  /* constant a [m] */
  double mRe;
  /* constant flat [-] */
  double mflat;
  /* constant omega [rad/s] */
  double momega;
  /* number of data sets */
  int mnum_data_sets;
  Datetime<nanoseconds> mtime_epoch;
  Datetime<nanoseconds> mfirst_epoch;
  Datetime<nanoseconds> mlast_epoch;

public:
  const char *agency() const noexcept { return charArena + agency_offset; }
  char *agency() noexcept { return charArena + agency_offset; }
  const char *satellite() const noexcept {
    return charArena + satellite_offset;
  }
  char *satellite() noexcept { return charArena + satellite_offset; }
  const char *sensor() const noexcept { return charArena + sensor_offset; }
  char *sensor() noexcept { return charArena + sensor_offset; }
  const char *process_level() const noexcept {
    return charArena + process_level_offset;
  }
  char *process_level() noexcept { return charArena + process_level_offset; }
  const char *pressure_type() const noexcept {
    return charArena + pressure_type_offset;
  }
  char *pressure_type() noexcept { return charArena + pressure_type_offset; }
  int file_type() const noexcept { return mfile_type; }
  int &file_type() noexcept { return mfile_type; }
  int file_format() const noexcept { return mfile_format; }
  int &file_format() noexcept { return mfile_format; }
  int num_header_records() const noexcept { return mnum_header_records; }
  int &num_header_records() noexcept { return mnum_header_records; }
  int num_data_records() const noexcept { return mnum_data_records; }
  int &num_data_records() noexcept { return mnum_data_records; }
  int max_degree() const noexcept { return mmax_degree; }
  int &max_degree() noexcept { return mmax_degree; }
  int coeff_errors() const noexcept { return mcoeff_errors; }
  int &coeff_errors() noexcept { return mcoeff_errors; }
  int coeff_normalized() const noexcept { return mcoeff_normalized; }
  int &coeff_normalized() noexcept { return mcoeff_normalized; }
  int num_data_sets() const noexcept { return mnum_data_sets; }
  int &num_data_sets() noexcept { return mnum_data_sets; }
  double GM() const noexcept { return mGM; }
  double &GM() noexcept { return mGM; }
  double Re() const noexcept { return mRe; }
  double &Re() noexcept { return mRe; }
  double flat() const noexcept { return mflat; }
  double &flat() noexcept { return mflat; }
  double omega() const noexcept { return momega; }
  double &omega() noexcept { return momega; }
  Datetime<nanoseconds> time_epoch() const noexcept { return mtime_epoch; }
  Datetime<nanoseconds> &time_epoch() noexcept { return mtime_epoch; }
  Datetime<nanoseconds> first_epoch() const noexcept { return mfirst_epoch; }
  Datetime<nanoseconds> &first_epoch() noexcept { return mfirst_epoch; }
  Datetime<nanoseconds> last_epoch() const noexcept { return mlast_epoch; }
  Datetime<nanoseconds> &last_epoch() noexcept { return mlast_epoch; }

  Aod1bIn(const char *fn) : mfn(fn) {
    std::memset(charArena, '\0', MaxArenaChars);
  }

  int read_header() noexcept;
}; /* class Aod1bIn */

} /* namespace dso */

#endif
