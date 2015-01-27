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
  double t = 0.07995893223819302e0;
  double l,lp,f,d,om;
  iers2010::fundarg (t,l,lp,f,d,om);
  printf ("\n\t|dl|  = %20.18f",fabs (2.291187512612069099e0-l));
  printf ("\n\t|dlp| = %20.18f",fabs (6.212931111003726414e0-lp));
  printf ("\n\t|df|  = %20.18f",fabs (3.658025792050572989e0-f));
  printf ("\n\t|dd|  = %20.18f",fabs (4.554139562402433228e0-d));
  printf ("\n\t|dom| = %20.18f",fabs (-0.5167379217231804489e0-om));
  
  // subroutine FCNNUT
  printf ("\nFunction FCNNUT");
  printf ("\n\tAbs. differences in microarcseconds:");
  double mjd = 54790e0;
  double x,y,dx,dy;
  iers2010::fcnnut (mjd,x,y,dx,dy);
  printf ("\n\t|dx|  = %20.18f",fabs (-176.8012290066270680e0-x));
  printf ("\n\t|dy|  = %20.18f",fabs (-93.51855308903756736e0-y));
  printf ("\n\t|ddx| = %20.18f",fabs (3.745573770491803067e0-dx));
  printf ("\n\t|ddy| = %20.18f",fabs (3.745573770491803067e0-dy));
  
  // subroutine PMSDNUT2
  printf ("\nFunction PMSDNUT2");
  printf ("\n\tAbs. differences in microarcseconds:");
  double pm[2];
  double rmjd = 54335e0;
  iers2010::pmsdnut2 (rmjd,pm);
  printf ("\n\t|pm(1)|  = %20.18f",fabs (24.83144238273364834e0-pm[0]));
  printf ("\n\t|pm(2)|  = %20.18f",fabs (-14.09240692041837661e0-pm[1]));
  
  
  // subroutine UTLIBR
  printf ("\nFunction UTLIBR (2 tests)");
  printf ("\n\tAbs. differences in mus and mus/day:");
  rmjd = 44239.1e0;
  double dut1, dlod;
  iers2010::utlibr (rmjd,dut1,dlod);
  printf ("\n\t->test A.");
  printf ("\n\t|ddut1|  = %20.18f",fabs (2.441143834386761746e0-dut1));
  printf ("\n\t|ddut1|  = %20.18f",fabs (-14.78971247349449492e0-dlod));
  rmjd = 55227.4e0;
  iers2010::utlibr (rmjd,dut1,dlod);
  printf ("\n\t->test B.");
  printf ("\n\t|ddut1|  = %20.18f",fabs (-2.655705844335680244e0-dut1));
  printf ("\n\t|ddut1|  = %20.18f",fabs (27.39445826599846967e0-dlod));
  
  // subroutine ARG2
  printf ("\nFunction ARG2");
  printf ("\n\tAbs. differences in radians:");
  double iangle[] = {2.849663065753787805e0,6.28318080000000023e0,
    4.926040134021299366e0,1.608450491115348768e0,2.375021572352622456e0,
    0.4746414933980958040e0,3.908159227647345801e0,2.551018561669245344e0,
    5.041990012540757959e0,4.206816878908014701e0,1.608463638294885811e0};
  int iyear =  2008;
  double day = 311.5e0;
  double angle[11];
  iers2010::arg2 (iyear,day,angle);
  for (int i = 0;i<11;i++) printf ("\n\t|dangle(%02i)|  = %20.18f",i+1, fabs(angle[i]-iangle[i]));
  
  // subroutine DEHANTTIDEINEL
  printf ("\nFunction DEHANTTIDEINEL");
  printf ("\n\tAbs. differences in meters :");
  double xsta[] = {4075578.385e0,931852.890e0,4801570.154e0};
  double xsun[] = {137859926952.015e0,54228127881.4350e0,23509422341.6960e0};
  double xmon[] = {-179996231.920342e0,-312468450.131567e0,-169288918.592160e0};
  double jc1    = 0.092799473402287e0,jc2 = .0e0;
  double xcor[3];
  iers2010::dehanttideinel (xsta,xsun,xmon,2009,4,13,.0e0,xcor);
  printf ("\n\tVersion 1 (UTC input):");
  printf ("\n\t|dx|  = %20.18f",0.7700420357108125891e-01-xcor[0]);
  printf ("\n\t|dy|  = %20.18f",0.6304056321824967613e-01-xcor[1]);
  printf ("\n\t|dz|  = %20.18f",0.5516568152597246810e-01-xcor[2]);
  iers2010::dehanttideinel (xsta,xsun,xmon,0.092799473402287e0,.0e0,xcor);
  printf ("\n\tVersion 2 (TT input):");
  printf ("\n\t|dx|  = %20.18f",0.7700420357108125891e-01-xcor[0]);
  printf ("\n\t|dy|  = %20.18f",0.6304056321824967613e-01-xcor[1]);
  printf ("\n\t|dz|  = %20.18f",0.5516568152597246810e-01-xcor[2]);
  
  // subroutine RG_ZONT2
  printf ("\nFunction RG_ZONT22");
  printf ("\n\tAbs. differences (see results) :");
  t = .07995893223819302e0;
  double dut, /*dlod,*/ domega;
  iers2010::rg_zont2 (t, dut, dlod, domega);
  printf ("\n\t|ddut|   = %20.18f (seconds)",7.983287678576557467E-002-dut);
  printf ("\n\t|ddlod|  = %20.18f (seconds/day)",5.035331113978199288E-005-dlod);
  printf ("\n\t|ddomega|= %20.18f (radians/second)", -4.249711616463017E-014-domega);
  
  // subroutine CNMTX
  printf ("\nFunction CNMTX (used by ORTHO_EOP)");
  printf ("\n\tAbs. differences (see results) :");
  double ht[12],h[12];
   *(h+0) = 15.35873641938967360e0;
   *(h+1) = 9.784941251812741214e0;
   *(h+2) = -5.520740128266865554e0;
   *(h+3) = 3.575314211234633888e0;
   *(h+4) = -13.93717453496387648e0;
   *(h+5) = -9.167400321705855504e0;
   *(h+6) = 5.532815475865292321e0;
   *(h+7) = 9.558741883500834646e0;
   *(h+8) = -10.22541212627272600e0;
   *(h+9)= 0.8367570529461261231e0;
   *(h+10)= 1.946355176475630611e0;
   *(h+11)= -13.55702062247304696e0;
  iers2010::cnmtx (54964.0e0,ht);
  for (int i=0;i<12;i++) printf ("\n\t|dh(%02i)|   = %20.18f",i+1,fabs(h[i]-ht[i]));

  // subroutine ORTHO_EOP
  printf ("\nFunction ORTHO_EOP");
  printf ("\n\tAbs. differences in microarcseconds :");
  iers2010::ortho_eop (47100e0,h);
  printf ("\n\t|deop(1)|   = %20.18f",-162.8386373279636530e0-h[0]);
  printf ("\n\t|deop(2)|   = %20.18f",117.7907525842668974e0-h[1]);
  printf ("\n\t|deop(3)|   = %20.18f",-23.39092370609808214e0-h[2]);

  printf ("\n");
  return 0;
}
