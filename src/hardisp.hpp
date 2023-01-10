#ifndef __IERS_1010_HARDISP__
#define __IERS_1010_HARDISP__

#include "datetime/dtcalendar.hpp"
#include "doodson.hpp"
#include <cmath>
#include <vector>
#include <array>

namespace iers2010 {

constexpr const int NTIN = 11;
constexpr const int MAX_LINE_SZ = 256;

struct BlqSiteInfo {
  // Cartwright-Tayler numbers of tides used in Scherneck lists:
  // M2, S2, N2, K2, K1, O1, P1, Q1, Mf, Mm, Ssa
  // translated to Doodson Numbers:
  static constexpr std::array<dso::DoodsonNumber,NTIN> Doodsons = {
      {2, 0, 0, 0, 0, 0},  {2, 2, -2, 0, 0, 0}, {2, -1, 0, 1, 0, 0},
      {2, 2, 0, 0, 0, 0},  {1, 1, 0, 0, 0, 0},  {1, -1, 0, 0, 0, 0},
      {1, 1, -2, 0, 0, 0}, {1, -2, 0, 1, 0, 0}, {0, 2, 0, 0, 0, 0},
      {0, 1, 0, -1, 0, 0}, {0, 0, 2, 0, 0, 0}
  };
  char site[32] = {'\0'};      ///< site name
  double amplitudes[NTIN * 3]; ///< amplitudes [m], radial, west, south
  double phases[NTIN * 3];     ///< corresponding phase values [rad]

  int parse_site_name(const char *line) noexcept;
};// BlqSiteInfo

/// @brief Read a BLQ file and parse site-dependent info (amplittudes and 
///        phases)
int parse_blq(const char *blqfn,
                        std::vector<BlqSiteInfo> &blqInfoVec,
                        const std::vector<const char *> *sites) noexcept;

class Hardisp {
private:
  double beta[6], beta_freq[6];
public:
  int operator()(const dso::TwoPartDate &tt_mjd,
               const dso::TwoPartDate &ut1_mjd) noexcept;
  void tdfrph(const dso::DoodsonNumber &d, double &phase,
             double &freq) const noexcept {
    phase = d.phase(beta);
    freq = d.frequency(beta_freq);
  }
}; // Hardisp

namespace hisp {

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
