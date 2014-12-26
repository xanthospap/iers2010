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
  
 
  printf ("\n");
  return 0;
}
