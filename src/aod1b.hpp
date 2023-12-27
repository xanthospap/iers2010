/** @file
 * Utilities to read/parse AOD1B files
 */

#ifndef __DSO_AOD1B_IPARSER_GPP__
#define __DSO_AOD1B_IPARSER_GPP__

#include "datetime/calendar.hpp"

namespace dso {

class Aod1bIn {
private:

public:
  const char *agency() const noexcept;
  char &*agency() noexcept;
  const char *satellite() const noexcept;
  char &*satellite() noexcept;
  const char *sensor() const noexcept;
  char &*sensor() noexcept;
  int file_type() const noexcept {return mfile_type;}
  int &file_type() noexcept {return mfile_type;}
  int file_format() const noexcept {return mfile_format;}
  int &file_format() noexcept {return mfile_format;}
  int num_header_records() const noexcept {return mnum_header_records;}
  int &num_header_records() noexcept {return mnum_header_records;}
  int num_data_records() const noexcept {return mnum_data_records;}
  int &num_data_records() noexcept {return mnum_data_records;}
  Datetime<nanoseconds> time_epoch() const noexcept {return mtime_epoch;}
  Datetime<nanoseconds> &time_epoch() noexcept {return mtime_epoch;}
  Datetime<nanoseconds> first_epoch() const noexcept {return mfirst_epoch;}
  Datetime<nanoseconds> &first_epoch() noexcept {return mfirst_epoch;}
  Datetime<nanoseconds> last_epoch() const noexcept {return mlast_epoch;}
  Datetime<nanoseconds> &last_epoch() noexcept {return mlast_epoch;}
}; /* class Aod1bIn */

}/* namespace dso */

#endif
