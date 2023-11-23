/** @file
 * Basic handling of icgem files (holding gravity potential models).
 * Information can be found on the International Centre for Global Earth
 * Models (ICGEM) website, http://icgem.gfz-potsdam.de/home
 * see Ince, E. S., Barthelmes, F., Reißland, S., Elger, K., Förste, C.,
 * Flechtner, F., Schuh, H. (2019): ICGEM – 15 years of successful collection
 * and distribution of global gravitational models, associated services and
 * future plans.- Earth System Science Data, 11, pp.647 - 674, DOI :
 * http:// doi.org/10.5194/essd-11-647-2019.
 *
 * ICGEM format: http://icgem.gfz-potsdam.de/ICGEM-Format-2023.pdf
 */

#ifndef __ICGEM_POTENTIAL_IO_HPP__
#define __ICGEM_POTENTIAL_IO_HPP__

#include "datetime/dtcalendar.hpp"
#include "stokes_coefficients.hpp"
#include <cstring>
#include <fstream>

namespace dso {

/** @brief A class to assist the reading/parsing of ICGEM gravity models.
 *
 * Download a gfc file from http://icgem.gfz-potsdam.de/tom_longtime and parse
 * it via this class. Note that the implementation is still incomplete and
 * not all models/parameters are read (aka only parameters of type 'gfc' can be
 * parsed but some models also have more coefficients).
 */
class Icgem {
public:
  typedef std::ifstream::pos_type pos_type;
  using Datetime = dso::datetime<dso::nanoseconds>;

  enum ErrorModel : char {No, Calibrated, Formal, CalibratedAndFormal};
  enum DataEntryType : char {gfc, gfct, trnd, asin, acos};
  enum IcgemVersion : char { v10, v20 };
  struct DataEntry {
    DataEntryType key;
    int degree; /** degree of current entry term */
    int order; /** order of current entry term */
    double C, S; /* Cnm and Snm harmonic terms */
    /* for ICGEM format version-1, t0 is t */
    Datetime t0, t1;
    double period; /* in [years] */
  }; /* struct DataEntry */

private:
  /** the filename */
  std::string _filename;
  /** start of data section (within file) */
  pos_type data_section_pos{0};
  /** the format of the file; default is v1.0, otherwise it must be specified 
   * within the file (header section). This is resilved at construction, when 
   * the header is parsed.
   */
  IcgemVersion _version{IcgemVersion::v10};
  /** the product type */
  char _product_type[64] = {'\0'};
  /** model name: name of the model (usually the respective filename without 
   * the extension “.gfc”)
   */
  char _modelname[128] = {'\0'};
  /** tide system: either
   *   "zero_tide",
   *   "tide_free",
   *   “mean_tide” or
   *   "unknown" (default)
   */
  char _tide_system[64] = "unknown";
  /** normalization: either "fully_normalized" (=default) or "unnormalized" */
  char _norm[64] = "fully_normalized";
  /** errors: either 
   *   "no", 
   *   "calibrated", 
   *   “formal” or 
   *   both "calibrated_and_formal" 
   * errors are included
   */
  ErrorModel _errors;
  /** gravitational constant times mass of the earth [m^3s^-2] */
  double _earth_gravity_constant{0e0};
  /** reference radius of the spherical harmonic development [m] */
  double _radius{0e0};
  /** maximum degree of the spherical harmonic development */
  int _max_degree{0};
  
  /** @brief Read the file header and assign basic information. 
   *
   * If the parameter \p quiet_skip is true, then the function will not emmit 
   * warnings for header lines that were not 'stored' in any of the member 
   * variables. If the \p quiet_skip is false, then any line within the header 
   * section that is not 'recognized' will emit a warning message to STDERR.
   * Warnings are only emitted for lines within the header section (i.e. 
   * between the begining of the file and 'end_of_head' part).
   */
  int parse_header(bool quiet_skip=true) noexcept;
  
  /** @brief Parse harmonic coefficients and store them in a StokesCoeffs 
   * instance.
   *
   * This function handles ICGEM file of format/version 2.0. It assumes that 
   * the file header has already been read (which is always true in case the 
   * instance has been successefuly created).
   * The function will read, parse and compute all Stokes coefficients Cnm, 
   * Snm up to degree \p up_to_degree and up to order \p up_to_order. 
   * Obviously, the \p up_to_degree must be less than or equal to the recorded 
   * 'max_degree' and greater than or equal to \p up_to_order.
   * The input epoch \p t is the epoch for which we want the Stokes 
   * Coefficients. It is needed to compute the time-varying part of the model, 
   * and/or validity intervals (parameters 'gfct', 'trnd', 'asin' and 'acos').
   * For each of the Stokes Coefficients, Cnm and Snm and for every degree and 
   * order pair (n,m), the coefficient are computed as:
   * C_nm(t) = C_nm(t0) + trnd*(t-t0) 
   *           + asin1*sin(2pi/p1 * (t-t0)) + acos1*cos(2pi/p1 * (t-t0))
   *           + asin2*sin(2pi/p2 * (t-t0)) + acos2*cos(2pi/p2 * (t-t0))
   * where, C_nm(t0) = gfc_nm + gfct_nm(t0),
   * and gfct, trnd, asin and acos are valid within the interval:
   * t0 <= t < t1
   * The same formula is used for the S_nm coefficients.
   * More infomation can be sought at the reference documents: 
   *  - http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf and 
   *  - http://icgem.gfz-potsdam.de/ICGEM-Format-2023.pdf
   *
   * @param[in] up_to_degree Max degree of S/C harmonic coeffients values to 
   *            read and store. The parameter must be <= coeffs.max_degree()
   *            and <= of the maximum degree recorded in the header section 
   *            (i.e. this.max_degree()).
   * @param[in] up_to_order Max order of S/C harmonic coeffients values to 
   *            read and store. The parameter must be <= coeffs.max_order() 
   *            and <= up_to_degree.
   * @param[in] t Epoch of computation (for the computed/stored Cnm and Snm 
   *            coefficients). Irrelevant if the model only contains static 
   *            terms, i.e. 'gfc'.
   * @param[out] coeffs An instance of type StokesCoeffs where the
   *            S/C harmonic coefficients are to be stored. If more space is 
   *            needed to store the coefficients, the instance will 
   *            automatically grow (i.e. no need to have the required capacity
   *            at input).
   */
  int parse_data_v2(int up_to_degree, int up_to_order, const Icgem::Datetime &t,
                    StokesCoeffs &coeffs) noexcept;
  
  /** @brief Parse harmonic coefficients and store them in a StokesCoeffs 
   * instance.
   *
   * This function handles ICGEM file of format/version 1.0. It assumes that 
   * the file header has already been read (which is always true in case the 
   * instance has been successefuly created).
   * The function will read, parse and compute all Stokes coefficients Cnm, 
   * Snm up to degree \p up_to_degree and up to order \p up_to_order. 
   * Obviously, the \p up_to_degree must be less than or equal to the recorded 
   * 'max_degree' and greater than or equal to \p up_to_order.
   * The input epoch \p t is the epoch for which we want the Stokes 
   * Coefficients. It is needed to compute the time-varying part of the model, 
   * and/or validity intervals (parameters 'gfct', 'trnd', 'asin' and 'acos').
   * For each of the Stokes Coefficients, Cnm and Snm and for every degree and 
   * order pair (n,m), the coefficient are computed as:
   * C_nm(t) = C_nm(t0) + trnd*(t-t0) 
   *           + asin1*sin(2pi/p1 * (t-t0)) + acos1*cos(2pi/p1 * (t-t0))
   *           + asin2*sin(2pi/p2 * (t-t0)) + acos2*cos(2pi/p2 * (t-t0))
   * where, C_nm(t0) = gfc_nm + gfct_nm(t0). 
   * In version 1.0 there are no validity intervals (i.e. t0 and t1 for each 
   * dynamic paramter) and the t0 values are given only for the 'gfct' terms.
   * Hence, 'asin', 'acos' and 'trnd' are supposed to hold for the same t0 
   * value as the 'gfct' term of the same (n,m) pair.
   *
   * The same formula is used for the S_nm coefficients.
   * More infomation can be sought at the reference documents: 
   *  - http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf and 
   *  - http://icgem.gfz-potsdam.de/ICGEM-Format-2023.pdf
   *
   * @param[in] up_to_degree Max degree of S/C harmonic coeffients values to 
   *            read and store. The parameter must be <= coeffs.max_degree()
   *            and <= of the maximum degree recorded in the header section 
   *            (i.e. this.max_degree()).
   * @param[in] up_to_order Max order of S/C harmonic coeffients values to 
   *            read and store. The parameter must be <= coeffs.max_order() 
   *            and <= up_to_degree.
   * @param[in] t Epoch of computation (for the computed/stored Cnm and Snm 
   *            coefficients). Irrelevant if the model only contains static 
   *            terms, i.e. 'gfc'.
   * @param[out] coeffs An instance of type StokesCoeffs where the
   *            S/C harmonic coefficients are to be stored. If more space is 
   *            needed to store the coefficients, the instance will 
   *            automatically grow (i.e. no need to have the required capacity
   *            at input).
   */
  int parse_data_v1(int up_to_degree, int up_to_order, const Icgem::Datetime &t,
                    StokesCoeffs &coeffs) noexcept;

public:

  /** Constructor from file name (of the icgem file).
   *
   * Note that the constructor will parse the file's header and assign member
   * variables (via a call to parse_header). In case this fails, then it will 
   * throw.
   */
  Icgem(const char *fn);

  /** Copy construction not allowed ! */
  Icgem(const Icgem &) noexcept = delete;

  /** Assignment operator not allowed */
  Icgem &operator=(const Icgem &) noexcept = delete;

  /** Destructor (no-op) */
  ~Icgem() noexcept {};
  
  /** get the filename */
  const std::string &filename() const noexcept {
    return _filename;
  }
  
  /** get the product type */
  const char *product_type() const noexcept {return _product_type;}
  
  /** get the model name */
  const char *model_name() const noexcept {return _modelname;}
  
  /** return the tide system (of the model) */
  const char *tide_system() const noexcept {return _tide_system;}
  
  /** get the supplied error 'model' */
  ErrorModel errors() const noexcept {return _errors;}

  /** version of the file */
  IcgemVersion version() const noexcept {return _version;}

  /** get max harmonics degree */
  int max_degree() const noexcept { return _max_degree; }

  /** get reference radius of the spherical harmonic development [m] */
  double radius() const noexcept { return _radius; }

  /** get gravitational constant times mass of the earth [m^3s^-2] */
  double gm() const noexcept { return _earth_gravity_constant; }

  /** return true if Stokes coefficients are normalized */
  bool is_normalized() const noexcept {
    return !std::strcmp(_norm, "fully_normalized");
  }
  
  /** @brief Parse harmonic coefficients and store them in a StokesCoeffs 
   * instance.
   *
   * This function handles ICGEM file of format/versions 1.0 and 2.0. It 
   * assumes that the file header has already been read (which is always true 
   * in case the instance has been successefuly created).
   * The function will read, parse and compute all Stokes coefficients Cnm, 
   * Snm up to degree \p up_to_degree and up to order \p up_to_order. 
   * Obviously, the \p up_to_degree must be less than or equal to the recorded 
   * 'max_degree' and greater than or equal to \p up_to_order.
   * The input epoch \p t is the epoch for which we want the Stokes 
   * Coefficients. It is needed to compute the time-varying part of the model, 
   * and/or validity intervals (parameters 'gfct', 'trnd', 'asin' and 'acos').
   * For each of the Stokes Coefficients, Cnm and Snm and for every degree and 
   * order pair (n,m), the coefficient are computed as:
   * C_nm(t) = C_nm(t0) + trnd*(t-t0) 
   *           + asin1*sin(2pi/p1 * (t-t0)) + acos1*cos(2pi/p1 * (t-t0))
   *           + asin2*sin(2pi/p2 * (t-t0)) + acos2*cos(2pi/p2 * (t-t0))
   * where, C_nm(t0) = gfc_nm + gfct_nm(t0).
   * For version 2.0 files, gfct, trnd, asin and acos should valid within the 
   * interval: t0 <= t < t1.
   * For version 1.0 files, there are no validity intervals (i.e. t0 and t1 
   * for each dynamic paramter) and the t0 values are given only for the 
   * 'gfct' terms. Hence, 'asin', 'acos' and 'trnd' are supposed to hold for 
   * the same t0 value as the 'gfct' term of the same (n,m) pair.
   *
   * The same formula is used for the S_nm coefficients.
   * More infomation can be sought at the reference documents: 
   *  - http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf and 
   *  - http://icgem.gfz-potsdam.de/ICGEM-Format-2023.pdf
   *
   * @param[in] up_to_degree Max degree of S/C harmonic coeffients values to 
   *            read and store. The parameter must be <= coeffs.max_degree()
   *            and <= of the maximum degree recorded in the header section 
   *            (i.e. this.max_degree()).
   *            If the parameter is not set, or if it is < 0, then the function 
   *            assumes that we are extracting coefficients up to the file's 
   *            max_degree (i.e. this->max_degree()).
   *            In case of ambiguity, you can always prompt the resulting 
   *            StokesCoeffs instance for the actual max degree.
   * @param[in] up_to_order Max order of S/C harmonic coeffients values to 
   *            read and store. The parameter must be <= coeffs.max_order() 
   *            and <= up_to_degree.
   *            If the parameter is not set, or if it is < 0, then the function 
   *            assumes that we are extracting coefficients up to the file's 
   *            max order.
   *            In case of ambiguity, you can always prompt the resulting 
   *            StokesCoeffs instance for the actual max order.
   * @param[in] t Epoch of computation (for the computed/stored Cnm and Snm 
   *            coefficients). Irrelevant if the model only contains static 
   *            terms, i.e. 'gfc'.
   * @param[out] coeffs An instance of type StokesCoeffs where the
   *            S/C harmonic coefficients are to be stored. If more space is 
   *            needed to store the coefficients, the instance will 
   *            automatically grow (i.e. no need to have the required capacity
   *            at input).
   * @return Anythin other than zero denotes an error.
   */
  int parse_data(int max_degree, int max_order, const Icgem::Datetime &t,
                 StokesCoeffs &coeffs) noexcept {
    if (max_degree < 0)
      max_degree = this->max_degree();
    if (max_order < 0)
      max_order = max_degree;
    if (version() == IcgemVersion::v10)
      return parse_data_v1(max_degree, max_order, t, coeffs);
    else /*if (version() == IcgemVersion::v20)*/
      return parse_data_v2(max_degree, max_order, t, coeffs);
  }

  int parse_data(const Icgem::Datetime &t, StokesCoeffs &coeffs) noexcept {
    return parse_data(-1, -1, t, coeffs);
  }

}; /* class Icgem */

  /** Parse & resolve an ICGEM v1.0 data line.
   *
   * The line can hold any of the 'gfc', 'gfct', 'trnd', 'asin' and 'acos'
   * parameters of the model. The resolved line is stored in an Icgem::DataEntry
   * instance. The ErrorModel of the file is needed to account for a correct
   * parsing of the fields.
   *
   * @param[in] line An ICGEM v1.0 data line
   * @param[in] er   The error model of the corresponding file
   * @param[out] entry The parameters of the line stored as an Icgem::DataEntry
   *                 instance. Note that for the v1.0 format, there are no
   *                 t0 and t1 values, but only a t value (for the 'gfct'
   *                 parameter). This is tored in the t0 value.
   * @return Anything other than 0 denotes an error.
   */
  int resolve_icgem_data_line_v1(const char *line, Icgem::ErrorModel er,
                                 Icgem::DataEntry &entry) noexcept;

  /** Parse & resolve an ICGEM v2.0 data line.
   *
   * The line can hold any of the 'gfc', 'gfct', 'trnd', 'asin' and 'acos'
   * parameters of the model. The resolved line is stored in an Icgem::DataEntry
   * instance. The ErrorModel of the file is needed to account for a correct
   * parsing of the fields.
   *
   * @param[in] line An ICGEM v2.0 data line
   * @param[in] er   The error model of the corresponding file
   * @param[out] entry The parameters of the line stored as an Icgem::DataEntry
   *                 instance.
   * @return Anything other than 0 denotes an error.
   */
  int resolve_icgem_data_line_v2(const char *line, Icgem::ErrorModel er,
                                 Icgem::DataEntry &entry) noexcept;

  /** Read and parse a date given in the format: yyyymmdd.hhmm
   *
   * This format is normally used within the data section of the ICGEM v2 format
   * (for validity intervals).
   * It is expected that the input string str starts with the first digit of
   * year (i.e. no whitespaces before the date).
   *
   * @param[in] line  A null-terminated C-string of type: yyyymmdd.hhmm
   * @param[out] t    The resolved datetime
   * @param[out] last The first non-resolved character (i.e. one character past
   *                  the yyyymmdd.hhmm part)
   * @return Anything other than 0 denotes an error
   */
  int resolve_icgem_data_date_v2(const char *line, Icgem::Datetime &t,
                                 const char *&last) noexcept;

  /** Read and parse a date given in the format: yyyymmdd
   *
   * This format is normally used within the data section of the ICGEM v1
   * format.
   * It is expected that the input string str starts with the first digit of
   * year (i.e. no whitespaces before the date).
   *
   * @param[in] line  A null-terminated C-string of type: yyyymmdd
   * @param[out] t    The resolved datetime
   * @param[out] last The first non-resolved character (i.e. one character past
   *                  the yyyymmdd part)
   * @return Anything other than 0 denotes an error
   */
  int resolve_icgem_data_date_v1(const char *line, Icgem::Datetime &t,
                                 const char *&last) noexcept;

} /* namespace dso */
#endif
