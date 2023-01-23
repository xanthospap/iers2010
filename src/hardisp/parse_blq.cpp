#include "hardisp.hpp"
#include "iers2010.hpp"
#include <charconv>
#include <fstream>
#include <cassert>

using dso::BlqSiteInfo;
constexpr const int NTIN = dso::BlqSiteInfo::NTIN;

namespace {
/// @brief Resolve and copy a site name from a corresponding BLQ line to
/// the buf buffer. 
/// Note that we are considering a max length of characters for a site name
/// of 31 chars
/// @param[in] buf A max of 32 characters (inluding \0) are copied to 
///                the buffer, hence its size should be 32
int parse_site_name(const char *line, char *buf) noexcept {
  // beacuse the site name can contain whitespace characters, start from the
  // end of the line and copy as many charactes are needed but < 32
  const int sz = std::strlen(line);
  const char *end = line + sz - 1;
  while ((end > line) && *end == ' ')
    --end;
  // we could also have whitespace characters at the start of the line!
  const char *start = line;
  while (*start && *start==' ') ++start;
  const int cpchars = std::min(31, (int)(end - start)+1);
  assert(cpchars>=0);
  // std::memcpy(site, start, sizeof(char) * cpchars);
  std::memcpy(buf, start, sizeof(char) * cpchars);
  buf[cpchars] = '\0';
  return cpchars;
}

bool site_in_vector(const char *line,
                    const std::vector<const char *> *sites) noexcept {
  if (!sites)
    return true;
  // copy name of site in buf
  char buf[32];
  parse_site_name(line, buf);
  // compare
  auto it =
      std::find_if(sites->begin(), sites->end(), [=](const char *site_name) {
        return !std::strcmp(buf, site_name);
      });
  return (it != sites->end());
}
} // unnamed namespace

int BlqSiteInfo::parse_site_name(const char *line) noexcept {
  return ::parse_site_name(line, site);
}

int dso::parse_blq(const char *blqfn,
                        std::vector<BlqSiteInfo> &blqInfoVec,
                        const std::vector<const char *> *sites) noexcept {
  // clear result vector
  blqInfoVec.clear();

  // quick return
  if (sites)
    if (!sites->size()) return 0;

  std::ifstream blq(blqfn);
  if (!blq.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening BLQ file %s (traceback: %s)\n",
            blqfn, __func__);
    return 1;
  }

  // lambda function: skip whitespace characters
  auto skipws = [](const char *str) noexcept -> const char * {
    while (*str && *str == ' ')
      ++str;
    return str;
  };

  char line[MAX_LINE_SZ];
  int error = 0;
  int all_sites_collected = 0;
  while (blq.getline(line, MAX_LINE_SZ) && !error && !all_sites_collected) {
    // skip any line(s) starting with '$' and read in station name
    do {
      blq.getline(line, MAX_LINE_SZ);
    } while (line[0] == '$' && blq && (!error));
    // note, we could have reached EOF
    if (blq.eof()) break;
    if (site_in_vector(line, sites)) {
      BlqSiteInfo siteinfo;
      siteinfo.parse_site_name(line);

      // go to the next non-comment line, and read three consecutive lines
      // these are the amplitude coefficients (in radial, west, south). NTIN
      // values per line. Units [m]
      double d11[NTIN];
      for (int k = 0; k < 3; k++) {
        // skip any comment lines
        do {
          blq.getline(line, MAX_LINE_SZ);
        } while (line[0] == '$' && blq && (!error));
        const char *cpos = line;
        const int sz = std::strlen(line);
        // collect 11 values from this line
        for (int i = 0; i < NTIN; i++) {
          auto cres = std::from_chars(skipws(cpos), line + sz, d11[i]);
          error += (cres.ec != std::errc{});
          cpos = cres.ptr;
        }
        std::memcpy(siteinfo.amplitudes + k * NTIN, d11,
                    sizeof(double) * NTIN);
      }
      if (error) {
        fprintf(stderr,
                "[ERROR] Failed collecting amplitudes from BLQ file %s "
                "(traceback: %s)\n",
                blqfn, __func__);
        return 1;
      }
      // go to the next non-comment line, and read three consecutive lines
      // these are the phase coefficients (in radial, west, south). NTIN
      // values per line. Units [deg] in file, transformed to [rad]
      for (int k = 0; k < 3; k++) {
        // skip any comment lines
        do {
          blq.getline(line, MAX_LINE_SZ);
        } while (line[0] == '$' && blq && (!error));
        const char *cpos = line;
        const int sz = std::strlen(line);
        // collect 11 values from this line
        for (int i = 0; i < NTIN; i++) {
          auto cres = std::from_chars(skipws(cpos), line + sz, d11[i]);
          error += (cres.ec != std::errc{});
          cpos = cres.ptr;
        }
        for (int i=0; i<11; i++) d11[i] = dso::deg2rad(d11[i]);
        std::memcpy(siteinfo.phases + k * NTIN, d11, sizeof(double) * NTIN);
      }
      if (error) {
        fprintf(stderr,
                "[ERROR] Failed collecting phases from BLQ file %s (traceback: "
                "%s)\n",
                blqfn, __func__);
        return 1;
      }

      // add new info in vector
      blqInfoVec.push_back(siteinfo);
      // have we collected info for all sites ?
      all_sites_collected = (sites)?(blqInfoVec.size()==sites->size()):0;
    }
    // done for current site !
  }

  // is everything ok?
  if (!sites) {
#ifdef DEBUG
    if (!blq.eof()) {
      fprintf(stderr,
              "[ERROR] Should have reached EOF but haven't! Failed parsing BLQ "
              "file %s (traceback: %s)\n",
              blqfn, __func__);
    }
#endif
    return blq.eof() ? 0 : 1;
  } else {
#ifdef DEBUG
    if (!all_sites_collected) {
      fprintf(stderr,
              "[WRNNG] Should have collected %d site, but collected %d from "
              "BLQ file %s (traceback: %s)\n",
              (int)sites->size(), (int)blqInfoVec.size(), blqfn, __func__);
    }
#endif
    return /*!all_sites_collected*/0;
  }
}
