#ifndef __BLQ_STREAM_HPP__
#define __BLQ_STREAM_HPP__

#include <fstream>

namespace iers2010
{
class BlqIn
{
public:
  /// Constructor from filesname
  explicit
  BlqIn(const char*);

  /// Copy not allowed
  BlqIn(const BlqIn&) = delete;

  /// Assignment not allowed
  BlqIn& operator=(const BlqIn&) = delete;

  /// Move constructor
  BlqIn(BlqIn&& b) = default;

  /// Move assignment operator
  BlqIn& operator=(BlqIn&& b) = default;

  int
  peak_next_station(std::string& sta);

  int
  skip_next_station();

  int
  read_next_station(std::string& sta, double tamp[3][11], double tph[3][11], 
    bool change_phase_sign=false);

  bool
  find_station(const std::string& station);
private:
  typedef std::ifstream::pos_type pos_type;

  int
  read_header();

  void
  goto_eoh();

  std::string    __filename;
  std::ifstream  __istream;
  pos_type       __eoheader;

};// class BlqIn

} // namespace iers2010

#endif
