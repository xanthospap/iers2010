/** @file
 * Utilities to read/parse AOD1B files
 */

#ifndef __DSO_AOD1B_IPARSER_GPP__
#define __DSO_AOD1B_IPARSER_GPP__

#include "datetime/calendar.hpp"
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace dso {

/** Different sets of Stokes coefficients included in AOD1B products */
enum class AOD1BCoefficientType { ATM, OCN, GLO, OBA };

/** A class to assist parsing of AOD1B product files.
 * Documentation can be found at:
 * ftp://isdcftp.gfz-potsdam.de/grace/DOCUMENTS/Level-1/
 * Reference document:
 * GRACE_AOD1B_Product_Description_Document_for_RL06.pdf
 *
 * This class is only meant to parse/store/interface-with the header part of
 * the product file. It also provides a list of function that can assist the
 * parsing of data and especially navigating through its header/data blocks.
 * It does not hold a stream to the file though (i.e. an ifstream instance)
 * but can act as an intermediate for more robust and sophisticated streams,
 * or more precisely stream iterators to the AOD1B file. See for example the
 * class Aod1bDataBlockIterator.
 */
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

public:
  /** A struct to represent a data block header for any coefficient type */
  struct Aod1bBlockHeader {
    int mset_nr;
    int mnum_lines;
    Datetime<nanoseconds> mepoch;
    AOD1BCoefficientType mtype;
  }; /* Aod1bBlockHeader */

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

  /** Given a data block hedaer line, this function will parse it into a
   * Aod1bBlockHeader instance
   */
  int parse_data_block_header(const char *line,
                              Aod1bBlockHeader &rec) const noexcept;

  /** Given a stream (fin) to the top of an AOD1B file, skip the header block */
  int skip_header(std::ifstream &fin) const noexcept;

  /** Given a stream (fin) to an AOD1B file, go to the next data block header.
   * Note that the stream should be placed at a position so that the next line
   * to be read is a data block header
   */
  int goto_next_block(std::ifstream &fin, Aod1bBlockHeader &rec) const noexcept;

  /** Read and skip a given number of lines within an AOD1B file (starting
   * from the given position within the file)
   */
  int skip_lines(std::ifstream &fin, int num_lines) const noexcept;

  /** Given a stream (fin) to an AOD1B file, go to the start of the next data
   * block of type \p type.
   * Note that the stream should be placed at a position so that the next line
   * to be read is a data block header
   */
  int goto_next_block(std::ifstream &fin, AOD1BCoefficientType type,
                      Aod1bBlockHeader &rec) const noexcept;

  template <AOD1BCoefficientType T> friend class Aod1bDataBlockIterator;

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

  /** Constructor from AOD1B filename.
   * The constructor will automatically call read_header on the instance,
   * and collect/assign all header info.
   */
  Aod1bIn(const char *fn);

  Aod1bIn(const Aod1bIn &other) noexcept;
  Aod1bIn(Aod1bIn &&other) noexcept;
  Aod1bIn &operator=(const Aod1bIn &other) noexcept;
  Aod1bIn &operator=(Aod1bIn &&other) noexcept;
  ~Aod1bIn() noexcept {};

  /** Read and parse an AOD1B header block, assigning info to the instance */
  int read_header() noexcept;

}; /* class Aod1bIn */

/** A comfort class to assist iterating through an AOD1B file.
 * Iterating is meant here in the sense of iterating through the data blocks
 * of the file, for a given coefficient type (i.e. ATM, OCN, GLO, etc).
 */
template <AOD1BCoefficientType T> class Aod1bDataBlockIterator {
  /* header infor and filename of the AOD1B file */
  const Aod1bIn *maod;
  /* a stream to the file */
  std::ifstream mfin;
  /* current block header */
  Aod1bIn::Aod1bBlockHeader mheader;

public:
  /** Constructor; we need a Aod1bIn instance (and not e.g. the AOD1B 
   * filename) so that we know that all relevant info are already parsed 
   * off from the file.
   * We will also try to get a stream to the file.
   */
  Aod1bDataBlockIterator(const Aod1bIn &aod)
      : maod(&aod), mfin(aod.mfn.c_str()){};

  /** get current header */
  const Aod1bIn::Aod1bBlockHeader &header() const noexcept { return mheader; }

  /** set the instance to the first data block in the AOD1B file */
  Aod1bDataBlockIterator &set_begin() {
    /* goto top of file and skip the header part */
    if (maod->skip_header(mfin)) {
      throw std::runtime_error(
          "[ERROR] Failed to skip header off from AOD1B file\n");
    }
    /* read untill we meet a data block header of given type */
    int j = maod->goto_next_block(mfin, T, mheader);
    /* hopefully nothing went wrong ... */
    if (j < 0 || j > 0) {
      throw std::runtime_error(
          "[ERROR] Failed to get begin data block from AOD1B file\n");
    }
    return *this;
  }

  /** skip the current data block (i.e. the one corresponding to this mheader */
  void skip() noexcept { maod->skip_lines(mfin, mheader.mnum_lines); }

  /** go/advance to the next data block */
  int advance() noexcept { return maod->goto_next_block(mfin, T, mheader); }
}; /* Aod1bDataBlockIterator<T> */

} /* namespace dso */

#endif
