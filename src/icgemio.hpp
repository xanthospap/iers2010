/** @file
 * Basic handling of icgem files (holding gravity potential models).
 * Information can be found on the International Centre for Global Earth
 * Models (ICGEM) website, http://icgem.gfz-potsdam.de/home
 * see Ince, E. S., Barthelmes, F., Reißland, S., Elger, K., Förste, C.,
 * Flechtner, F., Schuh, H. (2019): ICGEM – 15 years of successful collection
 * and distribution of global gravitational models, associated services and
 * future plans.- Earth System Science Data, 11, pp.647 - 674, DOI :
 * http:// doi.org/10.5194/essd-11-647-2019.
 */

#ifndef __ICGEM_POTENTIAL_IO_HPP__
#define __ICGEM_POTENTIAL_IO_HPP__

#include "datetime/tpdate.hpp"
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
  enum ErrorModel : char {No, Calibrated, Formal, CalibratedAndFormal};
  enum DataEntryType : char {gfc, gfct, trnd, asin, acos};
  struct DataEntry {
    DataEntryType key;
    int degree; /** degree of current entry term */
    int order; /** order of current entry term */
    double C, S; /* Cnm and Snm harmonic terms */
    /* for ICGEM format version-1, t0 is t */
    TwoPartDate t0, t1;
    double period; /* in [years] */
  }; /* struct DataEntry */

private:
  /** the filename */
  std::string _filename;
  /** start of data section (within file) */
  pos_type data_section_pos{0};
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
   * between the 'begin_of_head' and 'end_of_head' part.
   */
  int parse_header(bool quiet_skip=true) noexcept;

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

  /** @brief Parse harmonic coefficients up to degree l and order m.
   * 
   * Note that only data/values with a key value of 'gfc' are read. Some
   * models however, include parameters with more keys.
   * The function will also assign the model's constants (aka GM and Re) as
   * well as the normalization status of the coefficients.
   * 
   * @see http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf and 
   *      http://icgem.gfz-potsdam.de/ICGEM-Format-2023.pdf
   *
   * @param[in] l Max degree of S/C harmonic coeffients values to read and
   *            store.
   * @param[in] k Max order of S/C harmonic coeffients values to read and
   *            store (k <= l)
   * @param[in] t This is only used for Time Variable Gravity (TVG), where
   *            there are temporal coefficients.
   * @param[out] coeffs Pointer an instance of type HarmonicCoeffs where the
   *            S/C harmonic coefficients are to be stored. Note that this
   *            instance should have been allocated with enough space.
   */
  //int parse_data(int l, int k, const dso::TwoPartDate &t,
  //               StokesCoeffs *coeffs) noexcept;
}; /* class Icgem */

int resolve_icgem_data_line_v1(const char *line, Icgem::ErrorModel er,
                               Icgem::DataEntry &entry) noexcept;
int resolve_icgem_data_line_v2(const char *line, Icgem::ErrorModel er,
                               Icgem::DataEntry &entry) noexcept;
int resolve_icgem_data_date_v2(const char *line, dso::TwoPartDate &t,
                               const char *&last) noexcept;

} /* namespace dso */
#endif
