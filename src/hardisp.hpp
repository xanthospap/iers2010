#ifndef __IERS_1010_HARDISP__
#define __IERS_1010_HARDISP__

#include "datetime/dtcalendar.hpp"
#include <cmath>
#include <vector>

namespace iers2010 {

constexpr const int NTIN = 11;
constexpr const int MAX_LINE_SZ = 256;

struct BlqSiteInfo {
  char site[32] = {'\0'};      ///< site name
  double amplitudes[NTIN * 3]; ///< amplitudes (meter), radial, west, south
  double phases[NTIN * 3];     ///< corresponding phase values (degrees)

  int parse_site_name(const char *line) noexcept;
};// BlqSiteInfo

/// @brief Read a BLQ file and parse site-dependent info (amplittudes and 
///        phases)
int parse_blq(const char *blqfn,
                        std::vector<BlqSiteInfo> &blqInfoVec,
                        const std::vector<const char *> *sites) noexcept;

namespace hisp {

struct Constituent {
  dso::DoodsonNumber;
  double amplitude;
};// Constituent

// Parameters below set the buffer size for computing the tides recursively

// nl: the number of harmonics used in the prediction
constexpr int nl{600};

// nt: this must also be set in the subroutine admint
constexpr int nt{342};

// ntin: the number of harmonics read in
constexpr int ntin{11};

double eval(double, int, const double *, const double *, const double *);

int recurs(double *, int, const double *, int, const double *, double *);

int shells(double *, int *, int) noexcept;

int spline(int, const double *, const double *, double *, double *);

int tdfrph(const int *, dso::datetime<dso::seconds>, double &, double &);

int admint(const double *, const double *, dso::datetime<dso::seconds>,
           double *, double *, double *, int, int &);

int read_hardisp_args(double tamp[3][ntin], double tph[3][ntin],
                      const char *filename = nullptr);

int hardisp_impl(int, double, double tamp[3][ntin], double tph[3][ntin],
                 dso::datetime<dso::seconds> epoch, double *du, double *ds,
                 double *dw);

} // hisp

} // iers2010

#endif
