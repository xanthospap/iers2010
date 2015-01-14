#include <stdio.h>
#include "iers2010.hpp"

// COMPILATION: g++ -Wall -std=c++11 -L../lib/ -I../inc/ test.cpp -liers2010

int main ()
{
  printf ("\n===========================================================");
  printf ("\n TESTING IERS 2010 ROUTINES                                ");
  printf ("\n===========================================================");
  
  // subroutine FUNDARG
  printf ("\nFunction FUNDARG");
  printf ("\n\tAbs. differences in radians:");
  /*
   *  Test case:
   *     given input: T = 0.07995893223819302 Julian centuries since J2000
   *                  (MJD = 54465)
   *     expected output:  L = 2.291187512612069099 radians
   *                       LP = 6.212931111003726414 radians
   *                       F = 3.658025792050572989 radians
   *                       D = 4.554139562402433228 radians
   *                       OM = -0.5167379217231804489 radians
   */
  double t = 0.07995893223819302;
  double l,lp,f,d,om;
  iers2010::fundarg (t,l,lp,f,d,om);
  /*printf ("\n\t|dl|  = %15.12f",fabs (2.291187512612069099-l));
  printf ("\n\t|dlp| = %15.12f",fabs (6.212931111003726414-lp));
  printf ("\n\t|df|  = %15.12f",fabs (3.658025792050572989-f));
  printf ("\n\t|dd|  = %15.12f",fabs (4.554139562402433228-d));
  printf ("\n\t|dom| = %15.12f",fabs (-0.5167379217231804489-om));*/
  printf ("\n\t|dl|  = %20.15f",fabs (2.291187512612069099-l));
  printf ("\n\t|dlp| = %20.15f",fabs (6.212931111003726414-lp));
  printf ("\n\t|df|  = %20.15f",fabs (3.658025792050572989-f));
  printf ("\n\t|dd|  = %20.15f",fabs (4.554139562402433228-d));
  printf ("\n\t|dom| = %20.15f",fabs (-0.5167379217231804489-om));
  
  // subroutine FCNNUT
  printf ("\nFunction FCNNUT");
  printf ("\n\tAbs. differences in microarcseconds:");
  /*
   *  Test case: (NOT UPDATED FOR 2013 TABLE, TO BE DONE)
   *     given input: MJD = 54790D0   Modified Julian Date, TDB
   *                  
   *     expected output:  X = -176.8012290066270680D0 microarcseconds
   *                       Y = -93.51855308903756736D0 microarcseconds
   *                       dX = 3.745573770491803067D0 microarcseconds
   *                       dY = 3.745573770491803067D0 microarcseconds
   */
  double mjd = 54790e0;
  double x,y,dx,dy;
  iers2010::fcnnut (mjd,x,y,dx,dy);
  printf ("\n\t|dx|  = %15.12f",fabs (-176.8012290066270680e0-x));
  printf ("\n\t|dy|  = %15.12f",fabs (-93.51855308903756736e0-y));
  printf ("\n\t|ddx| = %15.12f",fabs (3.745573770491803067e0-dx));
  printf ("\n\t|ddy| = %15.12f",fabs (3.745573770491803067e0-dy));
  
  // subroutine PMSDNUT2
  printf ("\nFunction PMSDNUT2");
  printf ("\n\tAbs. differences in microarcseconds:");
  /*
   * Test case:
   *     given input: rmjd = 54335D0 ( August 23, 2007 ) 
   *
   *     expected output: (dx) pm(1)  = 24.83144238273364834D0 microarcseconds
   *                      (dy) pm(2) = -14.09240692041837661D0 microarcseconds
   */
  double pm[2];
  double rmjd = 54335e0;
  iers2010::pmsdnut2 (rmjd,pm);
  printf ("\n\t|pm(1)|  = %15.12f",fabs (24.83144238273364834e0-pm[0]));
  printf ("\n\t|pm(2)|  = %15.12f",fabs (-14.09240692041837661e0-pm[1]));
  
  
  // subroutine UTLIBR
  printf ("\nFunction UTLIBR (2 tests)");
  printf ("\n\tAbs. differences in mus and mus/day:");
  /*
   *     given input:  rmjd_a = 44239.1 ( January 1, 1980 2:24.00 )
   *                   rmjd_b = 55227.4 ( January 31, 2010 9:35.59 )
   *
   *     expected output: dUT1_a =   2.441143834386761746D0 mus;
   *                      dLOD_a = -14.78971247349449492D0 mus / day
   *                      dUT1_b = - 2.655705844335680244D0 mus;
   *                      dLOD_b =  27.39445826599846967D0 mus / day
   */
  rmjd = 44239.1e0;
  double dut1, dlod;
  iers2010::utlibr (rmjd,dut1,dlod);
  printf ("\n\t->test A.");
  printf ("\n\t|ddut1|  = %15.12f",fabs (2.441143834386761746e0-dut1));
  printf ("\n\t|ddut1|  = %15.12f",fabs (-14.78971247349449492e0-dlod));
  rmjd = 55227.4e0;
  iers2010::utlibr (rmjd,dut1,dlod);
  printf ("\n\t->test B.");
  printf ("\n\t|ddut1|  = %15.12f",fabs (-2.655705844335680244e0-dut1));
  printf ("\n\t|ddut1|  = %15.12f",fabs (27.39445826599846967e0-dlod));
  
  // subroutine ARG2
  printf ("\nFunction ARG2");
  printf ("\n\tAbs. differences in radians:");
  /*
   *  Test case:
   *     given input: IYEAR = 2008 
   *                  DAY = 311.5 (November 6 Noon)
   *     expected output: ANGLE(1)  = 2.849663065753787805D0  rad
   *                      ANGLE(2)  = 6.28318080000000023D0   rad
   *                      ANGLE(3)  = 4.926040134021299366D0  rad
   *                      ANGLE(4)  = 1.608450491115348768D0  rad
   *                      ANGLE(5)  = 2.375021572352622456D0  rad
   *                      ANGLE(6)  = 0.4746414933980958040D0 rad
   *                      ANGLE(7)  = 3.908159227647345801D0  rad
   *                      ANGLE(8)  = 2.551018561669245344D0  rad
   *                      ANGLE(9)  = 5.041990012540757959D0  rad 
   *                      ANGLE(10) = 4.206816878908014701D0  rad 
   *                      ANGLE(11) = 1.608463638294885811D0  rad 
   */
  double iangle[] = {2.849663065753787805e0,6.28318080000000023e0,
    4.926040134021299366e0,1.608450491115348768e0,2.375021572352622456e0,
    0.4746414933980958040e0,3.908159227647345801e0,2.551018561669245344e0,
    5.041990012540757959e0,4.206816878908014701e0,1.608463638294885811e0};
  int iyear =  2008;
  double day = 311.5e0;
  double angle[11];
  iers2010::arg2 (iyear,day,angle);
  for (int i = 0;i<11;i++) printf ("\n\t|dangle(%02i)|  = %15.12f",i+1, fabs (angle[i]-iangle[i]));
  
  // subroutine DEHANTTIDEINEL
  printf ("\nFunction DEHANTTIDEINEL");
  printf ("\n\tAbs. differences in meters :");
  /*
   * Test case:
   *     given input: XSTA(1) = 4075578.385D0 meters
   *                  XSTA(2) =  931852.890D0 meters
   *                  XSTA(3) = 4801570.154D0 meters   
   *                  XSUN(1) = 137859926952.015D0 meters
   *                  XSUN(2) = 54228127881.4350D0 meters
   *                  XSUN(3) = 23509422341.6960D0 meters
   *                  XMON(1) = -179996231.920342D0 meters
   *                  XMON(2) = -312468450.131567D0 meters
   *                  XMON(3) = -169288918.592160D0 meters
   *                  YR      = 2009
   *                  MONTH   = 4
   *                  DAY     = 13
   *                  FHR     = 0.00D0 hour
   *                  (T=0.092799473402287)
   *                  
   *     expected output:  DXTIDE(1) = 0.7700420357108125891D-01 meters
   *                       DXTIDE(2) = 0.6304056321824967613D-01 meters
   *                       DXTIDE(3) = 0.5516568152597246810D-01 meters
   */
  double xsta[] = {4075578.385e0,931852.890e0,4801570.154e0};
  double xsun[] = {137859926952.015e0,54228127881.4350e0,23509422341.6960e0};
  double xmon[] = {-179996231.920342e0,-312468450.131567e0,-169288918.592160e0};
  double jc1    = 0.092799473402287e0,jc2 = .0e0;
  double xcor[3];
  iers2010::dehanttideinel (xsta,xsun,xmon,2009,4,13,.0e0,xcor);
  printf ("\n\tVersion 1 (UTC input):");
  printf ("\n\t|dx|  = %15.12f",0.7700420357108125891e-01-xcor[0]);
  printf ("\n\t|dy|  = %15.12f",0.6304056321824967613e-01-xcor[1]);
  printf ("\n\t|dz|  = %15.12f",0.5516568152597246810e-01-xcor[2]);
  iers2010::dehanttideinel (xsta,xsun,xmon,0.092799473402287e0,.0e0,xcor);
  printf ("\n\tVersion 2 (TT input):");
  printf ("\n\t|dx|  = %15.12f",0.7700420357108125891e-01-xcor[0]);
  printf ("\n\t|dy|  = %15.12f",0.6304056321824967613e-01-xcor[1]);
  printf ("\n\t|dz|  = %15.12f",0.5516568152597246810e-01-xcor[2]);
  
  // subroutine RG_ZONT2
  printf ("\nFunction RG_ZONT22");
  printf ("\n\tAbs. differences (see results) :");
  /*
   *    Test case:
   *     given input: T = .07995893223819302 Julian centuries since J2000
   *                  (MJD = 54465)
   *     expected output: DUT    =  7.983287678576557467E-002 seconds
   *                      DLOD   =  5.035331113978199288E-005 seconds / day
   *                      DOMEGA = -4.249711616463017E-014 radians / second
   */
  t = .07995893223819302e0;
  double dut, /*dlod,*/ domega;
  iers2010::rg_zont2 (t, dut, dlod, domega);
  printf ("\n\t|ddut|   = %15.12f (seconds)",7.983287678576557467E-002-dut);
  printf ("\n\t|ddlod|  = %15.12f (seconds/day)",5.035331113978199288E-005-dlod);
  printf ("\n\t|ddomega|= %15.12f (radians/second)", -4.249711616463017E-014-domega);

  printf ("\n");
  return 0;
}
