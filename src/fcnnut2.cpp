#include <cassert>
#include <cmath>
#include <algorithm>
#include "iers2010.hpp"
#ifdef USE_EXTERNAL_CONSTS
#include "gencon.hpp"
#endif

// Block data of amplitudes for X (microas)
struct tbl_entry 
{
  double date, mxc, mxs, msx;
  double xc() const noexcept {return  mxc;}
  double xs() const noexcept {return  mxs;}
  double sx() const noexcept {return  msx;}    
  double yc() const noexcept {return  mxs;}
  double ys() const noexcept {return -mxc;}
  double sy() const noexcept {return  msx;}
};

const
std::vector<tbl_entry> table = {
  {45700.e0,     4.55e0,   -36.58e0,    19.72e0}, /// 1984.0
  {46066.e0,  -141.82e0,  -105.35e0,    11.12e0}, /// 1985.0
  {46431.e0,  -246.56e0,  -170.21e0,     9.47e0}, /// 1986.0
  {46796.e0,  -281.89e0,  -159.24e0,     8.65e0}, /// 1987.0
  {47161.e0,  -255.05e0,   -43.58e0,     8.11e0}, /// 1988.0
  {47527.e0,  -210.46e0,   -88.56e0,     7.31e0}, /// 1989.0
  {47892.e0,  -187.79e0,   -57.35e0,     6.41e0}, /// 1990.0
  {48257.e0,  -163.01e0,    26.26e0,     5.52e0}, /// 1991.0
  {48622.e0,  -145.53e0,    44.65e0,     4.80e0}, /// 1992.0
  {48988.e0,  -145.12e0,    51.49e0,     5.95e0}, /// 1993.0
  {49353.e0,  -109.93e0,    16.87e0,     9.45e0}, /// 1994.0
  {49718.e0,   -87.30e0,     5.36e0,     8.25e0}, /// 1995.0
  {50083.e0,   -90.61e0,     1.52e0,     7.67e0}, /// 1996.0
  {50449.e0,   -94.73e0,    35.35e0,     4.40e0}, /// 1997.0
  {50814.e0,   -67.52e0,    27.57e0,     3.40e0}, /// 1998.0
  {51179.e0,   -44.11e0,   -14.31e0,     3.45e0}, /// 1999.0
  {51544.e0,     5.21e0,   -74.87e0,     3.26e0}, /// 2000.0
  {51910.e0,    70.37e0,  -129.66e0,     2.86e0}, /// 2001.0
  {52275.e0,    86.47e0,  -127.84e0,     2.75e0}, /// 2002.0
  {52640.e0,   110.44e0,   -42.73e0,     2.59e0}, /// 2003.0
  {53005.e0,   114.78e0,    -0.13e0,     2.53e0}, /// 2004.0
  {53371.e0,   132.96e0,    -4.78e0,     2.72e0}, /// 2005.0
  {53736.e0,   157.36e0,    28.63e0,     2.19e0}, /// 2006.0
  {54101.e0,   160.40e0,    58.87e0,     1.87e0}, /// 2007.0
  {54466.e0,   156.76e0,   101.24e0,     1.74e0}, /// 2008.0
  {54832.e0,   142.99e0,   143.01e0,     1.89e0}, /// 2009.0
  {55197.e0,    33.70e0,   184.46e0,     1.95e0}, /// 2010.0
  {55562.e0,     0.76e0,   253.70e0,     1.14e0}, /// 2011.0
  {55927.e0,    25.47e0,   271.66e0,     1.07e0}, /// 2012.0
  {56293.e0,   113.42e0,   256.50e0,     1.86e0}  /// 2013.0
};
const int N = table.size();

/// @details This subroutine computes the effects of the free core nutation.  
///          Please note that the table is updated each year (see Note 4).
///          The parameter N needs to be incremented for each additional
///          year in the table.
/// 
/// @param[in] mjd   Modified Julian Date, TDB (Note 1)
/// @param[out]  x   CIP offset x component, in microas (Note 2)
/// @param[out]  y   CIP offset y component, in microas (Note 2)
/// @param[out] dx   Uncertainty of x component, in microas (Note 3)
/// @param[out] dy   Uncertainty of y component, in microas (Note 3) 
/// @return          An integer value, always 0.
///  
/// @warning Contains data table to be updated each year, 
///          see http://syrte.obspm.fr/~lambert/fcn/
/// 
/// @note 
///      -#  Though the Modified Julian Date (MJD) is strictly TDB, it is
///          usually more convenient to use TT, which makes no significant
///          difference.
///      -#  CIP is the Celestial Intermediate Pole. The expression
///          used is given in microarcseconds.
///      -#  The expression used is given in microarcseconds.
///      -#  The updated table is maintained at the website
///          http://syrte.obspm.fr/~lambert/fcn/.
/// 
/// @version 19.12.2013
/// 
/// @cite Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS 
///       Technical Note No. 36, BKG (2010); Chapter 5.5.5
int
iers2010::fcnnut2(double mjd, double& x, double& y, double& dx, double& dy) 
{

#ifdef USE_EXTERNAL_CONSTS
  constexpr double PI   (DPI);
#else
  constexpr double PI   (3.14159265358979323846e0);
#endif

  // Mean prediction error in microarcseconds per day
  constexpr double mpe    { 0.1325e0 };
  // FCN parameters period in days
  constexpr double per    { -430.21e0 };
  // Phase in radians
  double phi { (2e0*PI/per)*(mjd-51544.5e0) };


  // Prediction of the amplitude at the input date
  double axc {0e0},
         axs {0e0},
         ayc {0e0},
         ays {0e0};

  // First element in table for which date is larger than mjd
  auto it = std::upper_bound(table.begin(), table.end(), mjd,
    [](double d, const tbl_entry& t){return d<t.date;});

  if (it==table.begin()) {
    axc = it->xc();
    axs = it->xs();
    ayc = it->yc();
    ays = it->ys();
  } else if (it==table.end()) {
    axc = table[N-1].xc();
    axs = table[N-1].xs();
    ayc = table[N-1].yc();
    ays = table[N-1].ys();
  } else {
    --it;
    double t    = mjd - it->date;
    double dt   = (it+1)->date - it->date;
    double daxc = (it+1)->xc() - it->xc();
    double daxs = (it+1)->xs() - it->xs();
    double dayc = (it+1)->yc() - it->yc();
    double days = (it+1)->ys() - it->ys();
    axc  = it->xc() + (daxc/dt)*t;
    axs  = it->xs() + (daxs/dt)*t;
    ayc  = it->yc() + (dayc/dt)*t;
    ays  = it->ys() + (days/dt)*t;
  }

  // Computation of X and Y
  double cosphi {cos(phi)},
         sinphi {sin(phi)};
  x = axc*cosphi - axs*sinphi;  /// microas
  y = ayc*cosphi - ays*sinphi;  /// microas

  // Prediction of the uncertainty at the input date
  if (it==table.begin()) {
    axc = it->sx() + mpe*(it->date - mjd);
    axs = it->sx() + mpe*(it->date - mjd);
    ayc = it->sy() + mpe*(it->date - mjd);
    ays = it->sy() + mpe*(it->date - mjd);
  } else if (it==table.end()) {
    it = table.end()-1;
    axc = it->sx() + mpe*(mjd - it->date);
    axs = it->sx() + mpe*(mjd - it->date);
    ayc = it->sy() + mpe*(mjd - it->date);
    ays = it->sy() + mpe*(mjd - it->date);
  } else {
    double t    = mjd - it->date;
    double dt   = (it+1)->date - it->date;
    double daxc = (it+1)->sx() - it->sx();
    double daxs = (it+1)->sx() - it->sx();
    double dayc = (it+1)->sy() - it->sy();
    double days = (it+1)->sy() - it->sy();
    axc = std::abs(it->sx() + (daxc /dt)*t);
    axs = std::abs(it->sx() + (daxs /dt)*t);
    ayc = std::abs(it->sy() + (dayc /dt)*t);
    ays = std::abs(it->sy() + (days /dt)*t);
  }

  // Computation of the uncertainties
  dx = axc + axs; // microas
  dy = ayc + ays; // microas

  // Finished
  return 0;
}
