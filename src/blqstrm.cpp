#include <cstring>
#include <iostream>
#include <cerrno>
#include <cstdlib>
#include "blqstrm.hpp"

constexpr int MAX_HEADER_CHARS = 256;
constexpr int MAX_HEADER_LINES = 500;

/// Check if given line marks the end of header; for this to happen, the
/// line should:
/// 1. Start with '$$'
/// 2. Include the string 'END HEADER', procedded only by a number of whitespace 
///    or '$' characters
/// The function will return true only if the above conditions hold
bool
is_end_of_header(const char* line)
{
  int lensize = std::strlen(line);
  if (lensize<12) return false; // aka '$$END HEADER' = 12 chars
  if (line[0]!='$' || line[1]!='$') return false;
  const char* str = line;
  while (*str==' ' || *str=='$') ++str;
  return !std::strncmp(str, "END HEADER", 10);
}

iers2010::BlqIn::BlqIn(const char* filename)
  : __filename(filename),
    __istream(filename, std::ios_base::in),
    __eoheader(0)
{
  if (this->read_header()) {
    throw std::runtime_error("BlqIn::BlqIn() failed");
  }
}

/// 
int
iers2010::BlqIn::read_header()
{
  char line[MAX_HEADER_CHARS];

  // The stream should be open by now!
  if (!(this->__istream.is_open())) {
    std::cerr<<"\n[ERROR] BlqIn::read_header(): Stream is closed!";
    return 10;
  }

  // Go to the top of the file.
  __istream.seekg(0);

  // Read lines until end of header
  // ----------------------------------------------------
  int dummy_it = 0;
  __istream.getline(line, MAX_HEADER_CHARS);
  while (dummy_it<MAX_HEADER_LINES && !is_end_of_header(line)) {
    __istream.getline(line, MAX_HEADER_CHARS);
    ++dummy_it;
  }

  if (dummy_it>=MAX_HEADER_LINES) {
    std::cerr<<"\n[ERROR] BlqIn:read_header() -> Could not find 'END HEADER'.";
    return 3;
  }

  this->__eoheader = __istream.tellg();

  return 0;
}

/// Read next station from (an already open) BlqIn stream
//  47 $$
//  48   ACOR
//  49 $$ FES2004_PP ID: 2013-08-22 09:18:52
//  50 $$ Computed by OLMPP by H G Scherneck, Onsala Space Observatory, 2013
//  51 $$ ACOR,                      RADI TANG  lon/lat:  351.6011   43.3644   66.957
//  52   .03571 .01196 .00768 .00308 .00416 .00122 .00137 .00058 .00052 .00030 .00026
//  53   .00524 .00170 .00115 .00045 .00050 .00031 .00017 .00010 .00003 .00002 .00002
//  54   .00597 .00218 .00117 .00055 .00044 .00034 .00014 .00012 .00008 .00004 .00003
//  55    -88.1  -58.9 -107.3  -66.3  -77.9  171.7  -78.9   89.8    6.7    7.2    1.4
//  56     71.5   96.3   53.5   88.6   75.2  -18.7   74.1  -79.0  170.4 -168.8 -177.7
//  57    -57.0  -21.9  -77.9  -26.5  -46.2 -160.8  -46.9  148.0   -6.0    0.7   -0.4
int
iers2010::BlqIn::peak_next_station(std::string& sta)
{
  char line[MAX_HEADER_CHARS];
  const char* str;
  pos_type cur_pos = __istream.tellg();
  while (__istream.getline(line, MAX_HEADER_CHARS)) {
    if (line[0]!='$') {
      str = line;
      while (*str==' ' || *str=='$') ++str;
      if (std::strlen(str)>=4) {
        sta = std::string(str, 4);
        __istream.seekg(cur_pos);
        return 0;
      }
    }
  }

  if (__istream.eof()) return -1;
  return 1;
}

int
iers2010::BlqIn::skip_next_station()
{
  char line[MAX_HEADER_CHARS];

  // find station name ignoring comment lines
  while (__istream.getline(line, MAX_HEADER_CHARS)) {
    if (line[0]!='$') {
      const char* str = line;
      while (*str==' ' || *str=='$') ++str;
      if (std::strlen(str)>=4) break;
    }
  }
  // read comments if any
  while (__istream.getline(line, MAX_HEADER_CHARS)) {
    if (line[0]!='$') break;
  }
  // read values (exactly 6 rows)
  for (int i=0; i<6; i++) {
    if (!__istream.getline(line, MAX_HEADER_CHARS)) {
      std::cerr<<"\n[ERROR] Expected data line, found: \""<<line<<"\"";
      break;
    }
  }

  if (__istream.eof()) return -1;
  if (!__istream.good()) return 1;
  return 0;
}

/// @warning       Note that the function will change sign for phase, to be 
///                negative for lags
int
iers2010::BlqIn::read_next_station(std::string& sta, double tamp[3][11],
  double tph[3][11], bool change_phase_sign)
{
  char line[MAX_HEADER_CHARS];

  // find station name ignoring comment lines
  while (__istream.getline(line, MAX_HEADER_CHARS)) {
    if (line[0]!='$') {
      const char* str = line;
      while (*str==' ' || *str=='$') ++str;
      if (std::strlen(str)>=4) {
        sta = std::string(str, 4);
        break;
      }
    }
  }
  // read comments if any
  while (__istream.getline(line, MAX_HEADER_CHARS)) {
    if (line[0]!='$') break;
  }
  // read values (exactly 6 rows)
  char* end;
  const char *p;
  double dbl, *ar;
  for (int i=0; i<6; i++) {
    if  (i<3) ar = tamp[i];
    else      ar = tph[i-3];
    p = line;
    for (int j=0; j<11; j++) {
      dbl = std::strtod(p, &end);
      if (errno == ERANGE){
        std::cerr<<"\n[ERROR] Failed to resolve blq values for station: "<<sta<<" at line and column: "<<i+1<<"/"<<j+1;
        errno = 0;
        return 1;
      }
      ar[j] = dbl;
      p = end;
    }
    if (!__istream.getline(line, MAX_HEADER_CHARS)) {
      std::cerr<<"\n[ERROR] Expected data line, found: \""<<line<<"\"";
      break;
    }
  }

  if (change_phase_sign) 
    for (int i=0; i<3; i++) for (int j=0; j<11; j++) tph[i][j] = -tph[i][j];

  if (__istream.eof()) return -1;
  if (!__istream.good()) return 1;
  return 0;
}

void
iers2010::BlqIn::goto_eoh()
{
  __istream.clear();
  __istream.seekg(__eoheader);
}

bool
iers2010::BlqIn::find_station(const std::string& station)
{
  this->goto_eoh();
  std::string cur_sta;
  while (!peak_next_station(cur_sta)) {
    if (cur_sta==station) return true;
    else skip_next_station();
  }

  return false;
}
