/** @file
 * Utilities to read/parse AOD1B files
 */

#ifndef __DSO_AOD1B_IPARSER_GPP__
#define __DSO_AOD1B_IPARSER_GPP__

#include "datetime/calendar.hpp"
#include "stokes_coefficients.hpp"
#include "doodson.hpp"
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <array>

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
 *
 * Note that the epochs in AOD1B files are normally given in GPSTime, but the 
 * class holds dates in TT; the parsers take care of time scale 
 * transformations when the relevant functions are callled and all instance 
 * members (that hold datetimes) are stored in TT.
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
  /** A struct to represent a data block header for any coefficient type. */
  struct Aod1bBlockHeader {
    /** Data set number (i.e. within this AOD1B file, for a given coeff type */
    int mset_nr;
    /** Number of lines for this block */
    int mnum_lines;
    /** Epoch of block in TT */ 
    Datetime<nanoseconds> mepoch;
    /** Type of coefficients */
    AOD1BCoefficientType mtype;
  }; /* Aod1bBlockHeader */

private:
  /** filename */
  std::string mfn;
  /** buffer to store all character/string information in an AOD1B header */
  char charArena[MaxArenaChars];
  /** file type -- should be 999 */
  int mfile_type;
  /** file format: 0=binary 1=ascii */
  int mfile_format;
  /** number of header records */
  int mnum_header_records;
  /** number of data records */
  int mnum_data_records;
  /** maximum degree for harmonic coeffs */
  int mmax_degree;
  /** coefficient errors (yes/no) */
  int mcoeff_errors;
  /** coeff. normalized (yes/no) */
  int mcoeff_normalized;
  /** constant gm [m^3/s^2] */
  double mGM;
  /** constant a [m] */
  double mRe;
  /** constant flat [-] */
  double mflat;
  /** constant omega [rad/s] */
  double momega;
  /** number of data sets */
  int mnum_data_sets;
  /** Reference epoch in TT */
  Datetime<nanoseconds> mtime_epoch;
  /** Time of first data block in TT */
  Datetime<nanoseconds> mfirst_epoch;
  /** Time of last data block in TT */
  Datetime<nanoseconds> mlast_epoch;
  /** Only for TIDAL files */
  DoodsonConstituent mdoodson;
  
  /** Read and parse an AOD1B header block, assigning info to the instance */
  int read_header(std::ifstream &fin) noexcept;

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
   * from the given position within the file).
   */
  int skip_lines(std::ifstream &fin, int num_lines) const noexcept;

  /** Collect Stokes coefficients.
   * The function assumes that the stream passed in is placed just before 
   * the start of a (new) data block, i.e. the previous line read was the 
   * data block header.
   * It will continue reading lines and store coefficients for all 
   * (n,m)<=(max_degree,max_order). We do not assume that the the (n,m) 
   * coefficients are stored in some ordelry fashion, hence we must have an 
   * end condition, i.e. "when do we stop reading more lines ?". The 
   * criterion is given by the \p max_lines parameter. Hence, the function 
   * will continue to consume lines untill we reach max_lines.
   */
  int collect_coeffs(std::ifstream &fin, int max_degree, int max_order,
                     int max_lines, StokesCoeffs &cs) const noexcept;

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
  const std::string &fn() const noexcept {return mfn;}
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

private:
  /** Only used for very special cases; does not belong to class API */ 
  Aod1bIn() noexcept {};

public:
  Aod1bIn(const Aod1bIn &other) noexcept;
  Aod1bIn(Aod1bIn &&other) noexcept;
  Aod1bIn &operator=(const Aod1bIn &other) noexcept;
  Aod1bIn &operator=(Aod1bIn &&other) noexcept;
  ~Aod1bIn() noexcept {};

  /** Read and parse an AOD1B header block, assigning info to the instance */
  int read_header() noexcept;


  /** Read in and parse coefficients from a tidal AOD1B file, for a given
   * tidal constituent.
   * The passed in AOD1B file (as filename in \p fn), should be a TIDAL 
   * product, i.e. an AOD1B file containing coefficients for a given tidal
   * wave/constituent.
   *  
   * @warning This function is only meant to be used for TIDAL AOD1B files
   */
  static int get_tidal_wave_coeffs(const char *fn, StokesCoeffs &cCs,
                                   StokesCoeffs &sCs, Aod1bIn *aod1b,
                                   int max_degree = -1,
                                   int max_order = -1) noexcept;

}; /* class Aod1bIn */

/** A comfort class to assist iterating through an AOD1B file.
 * Iterating is meant here in the sense of iterating through the data blocks
 * of the file, for a given coefficient type (i.e. ATM, OCN, GLO, etc).
 *
 * Example usage:
 * // create an AodbIn instance from an AOD1B file
 * Aod1bIn aod(AODB1_FILE_NAME);
 * // a large enough Coefficients instance to hold data if needed
 * StokesCoeffs cs(aod.max_degree(), aod.max_degree(), ...); 
 * // create a data-block iterator, bind to the instance, for coefficients of 
 * // type 'ATM'
 * Aod1bDataBlockIterator<AOD1BCoefficientType::ATM> it(aod);
 * // set the iterator to the first data block
 * it.set_begin();
 * int error = 0;
 * while (!error) {
 *   if (SOME_CONDITION) {
 *     // collect all coefficients of the current block
 *     it.collect(cs);
 *   } else if (SOME_CONDITION) {
 *     // collect coefficients up to degree and order 120
 *     it.collect(cs, 120, 100);
 *   } else {
 *     // we don't need to store the data; skip it
 *     it.skip();
 *   }
 *   // advance iterator to the next data block
 *   j = it.advance();
 * }
 * // if j < 0 at the end of loop, we encountered EOF, i.e. we read through
 * // all the data blocks in the AOD1B file
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

  /** get a const reference to the underlying Aod1bIn instance */
  const Aod1bIn &aod1b() const {return *maod;}

  /** get current header */
  const Aod1bIn::Aod1bBlockHeader &header() const noexcept { return mheader; }

  /** mark this instance as EOF */
  void set_eof() noexcept {
    mheader.mepoch = Datetime<nanoseconds>::min();
  }

  /** check if this instance is marked as EOF */
  bool is_eof() const noexcept { 
    return (mheader.mepoch == Datetime<nanoseconds>::min() && mfin.eof());
  }

  /** Set the instance to the first data block in the AOD1B file.
   *
   * At success, the instance's mheader will hold the header of the first 
   * data block of the AOD1B file (of type T, e.g. "ATM").
   * The stream will be placed in a convinient position for a subsequent call 
   * to either skip() or collect() the data block.
   * In case of error, the function will throw.
   */
  Aod1bDataBlockIterator &set_begin() {
    /* first clear the stream (maybe it saw EOF) */
    mfin.clear();
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

  /** Skip the current data block (i.e. the one corresponding to this mheader */
  void skip() noexcept { maod->skip_lines(mfin, mheader.mnum_lines); }

  /** Collect the coefficients of the current data block.
   * 
   * Collect the coefficients of the current data block, up to maximum degree 
   * and order (max_degree, max_order). If these parameters are set to 
   * negative integers, all of the coefficients will be collected (i.e. max 
   * degree and order will be as reported in the AOD1B header).
   * The coefficients will be stored in the passed in StokesCoeffs instance, 
   * which must be large enough to hold them.
   * On success, 0 is returned. Any other integer denotes an error.
   */
  int collect(StokesCoeffs &cs, int max_degree=-1, int max_order=-1) noexcept {
    /* set max degree/order to collect */
    if (max_degree < 0)
      max_degree = maod->max_degree();
    if (max_order < 0)
      max_order = maod->max_degree();
    /* collect coefficients */
    if (maod->collect_coeffs(mfin, max_degree, max_order, mheader.mnum_lines,
                             cs)) {
      return 1;
    }
    return 0;
  }

  /** Go/advance to the next data block.
   *
   * Note that before calling this function, either skip() or collect should 
   * have been called, depending on what we want to do with data.
   * @return An integer denoting:
   *         -1 : EOF encountered
   *          0 : Success
   *          1 : Error
   * */
  int advance() noexcept { 
    int j = maod->goto_next_block(mfin, T, mheader);
    if (j<0) set_eof();
    return j;
  }
}; /* Aod1bDataBlockIterator<T> */

} /* namespace dso */

#endif
