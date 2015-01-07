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
  printf ("\n\t|dl|  = %15.12f",fabs (2.291187512612069099-l));
  printf ("\n\t|dlp| = %15.12f",fabs (6.212931111003726414-lp));
  printf ("\n\t|df|  = %15.12f",fabs (3.658025792050572989-f));
  printf ("\n\t|dd|  = %15.12f",fabs (4.554139562402433228-d));
  printf ("\n\t|dom| = %15.12f",fabs (-0.5167379217231804489-om));
  
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
 
  printf ("\n");
  return 0;
}
