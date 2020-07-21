#ifndef __IERS_1010_HARDISP__
#define __IERS_1010_HARDISP__

#include "ggdatetime/dtcalendar.hpp"
#include <cmath>

namespace iers2010 {

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

int tdfrph(const int *, ngpt::datetime<ngpt::seconds>, double &, double &);

int admint(const double *, const double *, ngpt::datetime<ngpt::seconds>,
           double *, double *, double *, int, int &);

int read_hardisp_args(double tamp[3][ntin], double tph[3][ntin],
                      const char *filename = nullptr);

int hardisp_impl(int, double, double tamp[3][ntin], double tph[3][ntin],
                 ngpt::datetime<ngpt::seconds> epoch);

} // namespace hisp

} // namespace iers2010

#endif
