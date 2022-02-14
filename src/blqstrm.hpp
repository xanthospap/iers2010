#ifndef __BLQ_IN_STREAM_HPP__
#define __BLQ_IN_STREAM_HPP__

#include <fstream>

namespace iers2010 {

/// @brief A class to hold input BLQ-format file
class BlqIn {
public:
  /// @brief Constructor from filename
  /// @param[in] filename  The filename of the BLQ file
  explicit BlqIn(const char *);

  /// Copy not allowed
  BlqIn(const BlqIn &) = delete;

  /// Assignment not allowed
  BlqIn &operator=(const BlqIn &) = delete;

  /// Move constructor
  BlqIn(BlqIn &&b) = default;

  /// Move assignment operator
  BlqIn &operator=(BlqIn &&b) = default;

  /// Destructor (does nothing)
  ~BlqIn() noexcept {};

  /// @brief Read next station from (an already open) BlqIn stream
  /// @param[out] sta  If found, the name of the next station BLQ record
  /// @return  An integer, denoting the following:
  ///          * -1  : EOF encountered before encountering next station
  ///          *  0  : all ok; station stored in sta and stream position is set
  ///                  to reading record of station sta
  ///          *  1  : error
  int peak_next_station(std::string &sta) noexcept;

  /// @brief read and skip next station's record
  /// @return An integer denoting:
  ///          * -1  : EOF encountered before reading next station
  ///          *  0  : all ok; station records read and skipped
  ///          *  1  : error
  int skip_next_station();

  /// @brief collect values fo the next station in the stream
  /// @param[out] sta The name of the station for which we collected the values
  /// @param[out] tamp The read in amplitudes
  /// @param[out] tph  The read in phases
  /// @param[in] change_phase_sign Set to true to change the sign of the phases
  ///            (aka every element in the tph array)
  /// @note Some algorithms (e.g. hardisp) use the phases with negative sign,
  ///       to be negative for lags
  int read_next_station(std::string &sta, double tamp[3][11], double tph[3][11],
                        bool change_phase_sign = false);

  /// @brief Find a station in the BLQ file stream (aka current instance)
  /// @param[in] station Station name to search for
  /// @return An boolean denoting:
  ///          *  true  : all ok; station found
  ///          *  false : station not found in file
  bool find_station(const std::string &station);

private:
  typedef std::ifstream::pos_type pos_type;

  /// @brief Read a BLQ-format file header.
  /// @return An integer; anything other than zero denotes an error:
  ///         *  0 : All ok
  ///         * 10 : Stream is closed
  ///         *  3 : Read MAX_HEADER_LINES before encountering the actual end of
  ///         header
  int read_header() noexcept;

  /// @brief Set the stream (get) position to end of header (__eoheader)
  void goto_eoh();

  std::string __filename;
  std::ifstream __istream;
  pos_type __eoheader;
}; // BlqIn

} // iers2010

#endif
