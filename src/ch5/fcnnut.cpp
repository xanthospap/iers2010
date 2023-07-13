#include "iers2010.hpp"
#include "iersc.hpp"
#include <array>
#include <cmath>

/* Block data of amplitudes for X (microas). Note that according to the
 * IERS2010 standards (see 5.5.5), the following relationships hold:
 * Ys = -Xc and
 * Yc = Xs
 * where XC and XS are the amplitudes of the cosine and sine terms,
 * respectively (here denoted by mxc, mxs).
 * Date is provided as an MJD numeric value.
 */
struct tbl_entry {
  dso::TwoPartDate date;
  double mxc, mxs, msx;
  constexpr double xc() const noexcept { return mxc; }
  constexpr double ys() const noexcept { return -mxc; }
  constexpr double xs() const noexcept { return mxs; }
  constexpr double yc() const noexcept { return mxs; }
  constexpr double sx() const noexcept { return msx; }
  constexpr double sy() const noexcept { return msx; }
};

/* Block data of amplitudes for X (microas), Table 5.2 from IERS 2010 */
constexpr const std::array<tbl_entry, 30> table = {{
    {dso::TwoPartDate(45700e0, 0e0), 4.55e0, -36.58e0, 19.72e0},     /// 1984.0
    {dso::TwoPartDate(46066e0, 0e0), -141.82e0, -105.35e0, 11.12e0}, /// 1985.0
    {dso::TwoPartDate(46431e0, 0e0), -246.56e0, -170.21e0, 9.47e0},  /// 1986.0
    {dso::TwoPartDate(46796e0, 0e0), -281.89e0, -159.24e0, 8.65e0},  /// 1987.0
    {dso::TwoPartDate(47161e0, 0e0), -255.05e0, -43.58e0, 8.11e0},   /// 1988.0
    {dso::TwoPartDate(47527e0, 0e0), -210.46e0, -88.56e0, 7.31e0},   /// 1989.0
    {dso::TwoPartDate(47892e0, 0e0), -187.79e0, -57.35e0, 6.41e0},   /// 1990.0
    {dso::TwoPartDate(48257e0, 0e0), -163.01e0, 26.26e0, 5.52e0},    /// 1991.0
    {dso::TwoPartDate(48622e0, 0e0), -145.53e0, 44.65e0, 4.80e0},    /// 1992.0
    {dso::TwoPartDate(48988e0, 0e0), -145.12e0, 51.49e0, 5.95e0},    /// 1993.0
    {dso::TwoPartDate(49353e0, 0e0), -109.93e0, 16.87e0, 9.45e0},    /// 1994.0
    {dso::TwoPartDate(49718e0, 0e0), -87.30e0, 5.36e0, 8.25e0},      /// 1995.0
    {dso::TwoPartDate(50083e0, 0e0), -90.61e0, 1.52e0, 7.67e0},      /// 1996.0
    {dso::TwoPartDate(50449e0, 0e0), -94.73e0, 35.35e0, 4.40e0},     /// 1997.0
    {dso::TwoPartDate(50814e0, 0e0), -67.52e0, 27.57e0, 3.40e0},     /// 1998.0
    {dso::TwoPartDate(51179e0, 0e0), -44.11e0, -14.31e0, 3.45e0},    /// 1999.0
    {dso::TwoPartDate(51544e0, 0e0), 5.21e0, -74.87e0, 3.26e0},      /// 2000.0
    {dso::TwoPartDate(51910e0, 0e0), 70.37e0, -129.66e0, 2.86e0},    /// 2001.0
    {dso::TwoPartDate(52275e0, 0e0), 86.47e0, -127.84e0, 2.75e0},    /// 2002.0
    {dso::TwoPartDate(52640e0, 0e0), 110.44e0, -42.73e0, 2.59e0},    /// 2003.0
    {dso::TwoPartDate(53005e0, 0e0), 114.78e0, -0.13e0, 2.53e0},     /// 2004.0
    {dso::TwoPartDate(53371e0, 0e0), 132.96e0, -4.78e0, 2.72e0},     /// 2005.0
    {dso::TwoPartDate(53736e0, 0e0), 157.36e0, 28.63e0, 2.19e0},     /// 2006.0
    {dso::TwoPartDate(54101e0, 0e0), 160.40e0, 58.87e0, 1.87e0},     /// 2007.0
    {dso::TwoPartDate(54466e0, 0e0), 156.76e0, 101.24e0, 1.74e0},    /// 2008.0
    {dso::TwoPartDate(54832e0, 0e0), 142.99e0, 143.01e0, 1.89e0},    /// 2009.0
    {dso::TwoPartDate(55197e0, 0e0), 33.70e0, 184.46e0, 1.95e0},     /// 2010.0
    {dso::TwoPartDate(55562e0, 0e0), 0.76e0, 253.70e0, 1.14e0},      /// 2011.0
    {dso::TwoPartDate(55927e0, 0e0), 25.47e0, 271.66e0, 1.07e0},     /// 2012.0
    {dso::TwoPartDate(56293e0, 0e0), 113.42e0, 256.50e0, 1.86e0}     /// 2013.0
}};
constexpr const int N = table.size();

int iers2010::fcnnut(const dso::TwoPartDate &mjd_tt, double &x, double &y,
                     double &dx, double &dy) noexcept {

  /* Mean prediction error in microarcseconds per day */
  constexpr const double mpe = 0.1325e0;

  /* FCN parameters period in days */
  constexpr const double per = -430.21e0;
  
  /* t is given in days since J2000.0 */
  const double tdays = mjd_tt.diff<dso::DateTimeDifferenceType::FractionalDays>(
      dso::TwoPartDate(dso::j2000_mjd, 0e0));
  
  /* Phase in radians */
  const double phi = (2e0 * iers2010::DPI / per) * tdays;

  /* Prediction of the amplitude at the input date */
  double axc(0e0), axs(0e0), ayc(0e0), ays(0e0);
  double saxc(0e0), saxs(0e0), sayc(0e0), says(0e0);

  /* search for suitable interval in reverse order */
  if (mjd_tt >= table[N - 1].date) {
    /* mjd >= last date in table  */
    axc = table[N - 1].xc();
    axs = table[N - 1].xs();
    ayc = table[N - 1].yc();
    ays = table[N - 1].ys();
    const double dt =
        mjd_tt.diff<dso::DateTimeDifferenceType::FractionalDays>(
            table[N - 1].date);
    saxc = table[N - 1].sx() + mpe * dt;
    saxs = saxc;
    sayc = table[N - 1].sy() + mpe * dt;
    sayc = sayc;
  } else {
    auto upper = table.end() - 1;
    /* search in reverse order */
    while (upper > table.begin()) {
      auto lower = upper - 1;
      if (mjd_tt >= lower->date && mjd_tt < upper->date) {
        const double t =
            mjd_tt.diff<dso::DateTimeDifferenceType::FractionalDays>(
                lower->date);
        const double dt =
            upper->date.diff<dso::DateTimeDifferenceType::FractionalDays>(
                lower->date);
        const double daxc = upper->xc() - lower->xc();
        const double daxs = upper->xs() - lower->xs();
        const double dayc = upper->yc() - lower->yc();
        const double days = upper->ys() - lower->ys();
        /* coefficients for X and Y interpolation */
        axc = lower->xc() + (daxc / dt) * t;
        axs = lower->xs() + (daxs / dt) * t;
        ayc = lower->yc() + (dayc / dt) * t;
        ays = lower->ys() + (days / dt) * t;
        /* coefficients for uncertainties */
        const double sdaxc = upper->sx() - lower->sx();
        const double sdaxs = sdaxc;
        const double sdayc = upper->sy() - lower->sy();
        const double sdays = sdayc;
        saxc = std::abs(lower->sx() + (sdaxc / dt) * t);
        saxs = std::abs(lower->sx() + (sdaxs / dt) * t);
        sayc = std::abs(lower->sy() + (sdayc / dt) * t);
        says = std::abs(lower->sy() + (sdays / dt) * t);
        /* break loop */
        break;
      }
      --upper;
    }
    if (upper == table.begin()) {
      /* given date before first table entry */
#ifdef DEBUG
      if (mjd_tt >= table.begin()->date) {
        fprintf(stderr, "ERROR %.20e >= %.20e\n", mjd_tt.as_mjd(),
                table.begin()->date.as_mjd());
      }
      assert(mjd_tt < table.begin()->date);
#endif
      axc = table.begin()->xc();
      axs = table.begin()->xs();
      ayc = table.begin()->yc();
      ays = table.begin()->ys();
      const double dt =
          table.begin()->date.diff<dso::DateTimeDifferenceType::FractionalDays>(
              mjd_tt);
      saxc = table.begin()->sx() + mpe * dt;
      saxs = saxc;
      sayc = table.begin()->sy() + mpe * dt;
      sayc = sayc;
    }
  }

  /* Computation of X and Y */
  const double cosphi = std::cos(phi);
  const double sinphi = std::sin(phi);
  x = axc * cosphi - axs * sinphi; /* μas */
  y = ayc * cosphi - ays * sinphi; /* μas */

  /* Computation of the uncertainties */
  dx = saxc + saxs; /* μas */
  dy = sayc + says; /* μas */

  /* Finished */
  return 0;
}
