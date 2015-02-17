#include <stdio.h>
#include <string.h>
#include "iers2010.hpp"

#define pi  3.1415926535e0

// COMPILATION: g++ -Wall -std=c++11 -L../lib/ -I../inc/ test.cpp -liers10++

int main ()
{
  int status = 0;
  double dlat[] = { 48.20e0,  89.20e0, -89.20e0};
  double dlon[] = { 16.37e0,  16.37e0,  16.37e0};
  double hell[] = {156.00e0, 156.00e0, 156.00e0};
  double dmjd = 56141.e0;
  int nstat = sizeof(dlat) / sizeof(double);
  int it = 0;
  
  double mp[nstat], mt[nstat], mdt[nstat], me[nstat], 
  mah[nstat], maw[nstat], mundu[nstat];
  
  for (int i=0;i<nstat;i++) {
    dlat[i] *= pi/180.e0;
    dlon[i] *= pi/180.e0;
  }
  
  status = iers2010::gpt2 (
    dmjd,dlat,dlon,hell,nstat,mp,mt,mdt,me,mah,maw,mundu,it,
    "/home/xanthos/Software/iers2010/src/gpt2_5.grd");
    
  if (status) {
    printf ("\nERROR! gpt2 could not run. Error code: %01i",status);
    status = 1;
   } else {
     for (int i=0;i<nstat;i++) {
       printf ("\nResults for station nr %1i", i+1);
        printf ("\n\t|dp|      = %10.3f hPa",mp[i]);
       printf ("\n\t|dT|      = %10.3f Celsius",mt[i]);
       printf ("\n\t|ddT|     = %10.3f deg/km",mdt[i]);
       printf ("\n\t|de|      = %10.3f hPa",me[i]);
       printf ("\n\t|dah|     = %15.8f (unitless)",mah[i]);
       printf ("\n\t|daw|     = %15.8f (unitless)",maw[i]);
       printf ("\n\t|dundu|   = %10.3f meters\n",mundu[i]);
     }
    }

  char ch_status[6];
  if (status)
    strcpy (ch_status,"ERROR");
  else
    strcpy (ch_status,"OK");
  printf ("\nStatus : %s", ch_status);

  printf ("\n");
  return 0;
}