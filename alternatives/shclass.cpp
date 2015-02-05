#include "shclass.hpp"
#include <algorithm>
#include <stdio.h>
//#include <exception>

iers2010::SphericalHarmonics::SphericalHarmonics ()
:  mdegree(9), morder(9), mx(0e0), my(0e0), mz(0e0)
{
  //printf ("\n\tINITIALIZED");
  mv = new double* [mdegree+1];
  mw = new double* [mdegree+1];
  for (int i=0;i<=mdegree;i++) {
    mv[i] = new double [i+1];
    mw[i] = new double [i+1];
   }
}

iers2010::SphericalHarmonics::SphericalHarmonics (const double& x, const double& y, const double& z)
:  mdegree(9), morder(9), mx(x), my(y), mz(z)
  {
    ///printf ("\n\tINITIALIZED");
    mv = new double* [mdegree+1];
    mw = new double* [mdegree+1];
    for (int i=0;i<=mdegree;i++) {
      mv[i] = new double [i+1];
      mw[i] = new double [i+1];
    }
    
    this->populate ();
  }
  
iers2010::SphericalHarmonics::~SphericalHarmonics ()
{
  //printf ("\n\tDYING");
  for (int i=0;i<=mdegree;i++) {
    delete [] mv[i];
    delete [] mw[i];
   }
   delete [] mv;
   delete [] mw;
}

iers2010::SphericalHarmonics::SphericalHarmonics (const SphericalHarmonics& s)
:  mdegree(9), morder(9), mx(s.mx), my(s.my), mz(s.mz)
{
  //printf ("\n\tINITIALIZED");
  mv = new double* [mdegree+1];
  mw = new double* [mdegree+1];
  for (int i=0;i<=mdegree;i++) {
    mv[i] = new double [i+1];
    std::copy (s.mv[i], s.mv[i]+i+1, mv[i]);
    mw[i] = new double [i+1];
    std::copy (s.mw[i], s.mw[i]+i+1, mw[i]);
   }
}

iers2010::SphericalHarmonics& iers2010::SphericalHarmonics::operator = (const SphericalHarmonics& s)
{  //printf ("\n\tINITIALIZED");
  if (this != &s) {
    mx = s.mx;
    my = s.my;
    mz = s.mz;
    this->populate ();
  }
  return *this;
}

void iers2010::SphericalHarmonics::populate ()
{
  // Legendre polynomials
  mv[0][0] = 1e0;
  mw[0][0] = 0e0;
  mv[1][0] = mz * mv[0][0];
  mw[1][0] = 0e0;
  
  for (int n=1;n<mdegree;n++) {
    int N ( n + 1 );
    mv[n+1][0] = ( (2*N-1) * mz * mv[n][0] - (N-1) * mv[n-1][0] ) / (double) N;
    mw[n+1][0] = 0e0;
   }
   
  for (int m=0;m<morder;m++) {
    int M ( m + 1 );
    mv[m+1][m+1] = (double) (2*M-1) * ( mx*mv[m][m] - my*mw[m][m] );
    mw[m+1][m+1] = (double) (2*M-1) * ( mx*mw[m][m] + my*mv[m][m] );
    if (m<morder-1) {
      mv[m+2][m+1] = (2*M+1) * mz* mv[m+1][m+1];
      mw[m+2][m+1] = (2*M+1) * mz* mw[m+1][m+1];
     }
    int N = M + 2;
    for (int n=m+2;n<mdegree;n++) {
      mv[n+1][m+1] = ( (2*N-1)*mz*mv[n][m+1] - (N+M-1)*mv[n-1][m+1] ) 
      / (double) (N-M);
      mw[n+1][m+1] = ( (2*N-1)*mz*mw[n][m+1] - (N+M-1)*mw[n-1][m+1] ) 
      / (double) (N-M);
      N++;
     }
   }
   
   return;
}