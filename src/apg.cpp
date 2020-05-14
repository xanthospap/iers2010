#include "iers2010.hpp"

/// @details  This subroutine determines the asymmetric delay d in meters caused
///           by gradients.  The north and east gradients are also provided.
///           They are based on Spherical Harmonics up to degree and order 9.
///           If the north and east gradients are used, they should be used
///           with the gradient model by Chen and Herring (1997). See Reference 1
///           and last lines of this subroutine.
///           This function is a translation/wrapper for the fortran APG
///           subroutine, found here : 
///           http://maia.usno.navy.mil/conv2010/software.html
/// 
/// @param[in]  dlat  Latitude given in radians (North Latitude)
/// @param[in]  dlon  Longitude given in radians (East Longitude)
/// @param[in]  az    Azimuth from north in radians
/// @param[in]  el    Elevation angle in radians
/// @param[out] d     Delay in meters
/// @param[out] grn   North gradient in mm
/// @param[out] gre   East gradient in mm
/// @return           An integer 
/// 
/// @note
///    -# This a priori model cannot replace the (additional) estimation of
///       gradient parameters, if observations at elevation angles below
///       15 degrees are analyzed.
///    -# Status: Class 1 model
///
/// @verbatim
/// Test case:
///      Kashima 11 Station information retrieved at:
///      ftp://ivscc.gsfc.nasa.gov/pub/config/ns/kashim11.config.txt
/// 
///      given input: DLAT = 0.6274877539940092D0 radians (KASHIMA 11, Japan)
///                   DLON = 2.454994088489240D0 radians
///                   AZ   = 0.2617993877991494D0 radians
///                   EL   = 0.8726646259971648D0 radians
/// 
///      expected output: D   = -0.9677190006296187757D-4 meters 
///                       GRN = -0.1042668498001996791D0 mm
///                       GRE = 0.4662515377110782594D-1 mm
/// @endverbatim
/// 
/// @version 29.09.2010
/// 
/// @cite iers2010
///     Chen, G. and Herring, T. A., 1997, ``Effects of atmospheric azimuthal
///     asymmetry on the analysis of space geodetic data,"
///     J. Geophys. Res., 102(B9), pp. 20,489--20,502, doi: 10.1029/97JB01739.
/// 
int
iers2010::apg(double dlat, double dlon, double az, double el, 
  double& d, double& grn, double& gre)
{
  // degree n and order m
  constexpr int nmax = 9;
  constexpr int mmax = 9;
    
  /*static*/ constexpr double a_n[] = {
     2.8959e-02,-4.6440e-01,-8.6531e-03, 1.1836e-01,-2.4168e-02,
    -6.9072e-05, 2.6783e-01,-1.1697e-03,-2.3396e-03,-1.6206e-03,
    -7.4883e-02, 1.3583e-02, 1.7750e-03, 3.2496e-04, 8.8051e-05,
     9.6532e-02, 1.3192e-02, 5.5250e-04, 4.0507e-04,-5.4758e-06,
     9.4260e-06,-1.0872e-01, 5.7551e-03, 5.3986e-05,-2.3753e-04,
    -3.8241e-05, 1.7377e-06,-4.4135e-08, 2.1863e-01, 2.0228e-02,
    -2.0127e-04,-3.3669e-04, 8.7575e-06, 7.0461e-07,-4.0001e-08,
    -4.5911e-08,-3.1945e-03,-5.1369e-03, 3.0684e-04, 2.4459e-05,
     7.6575e-06,-5.5319e-07, 3.5133e-08, 1.1074e-08, 3.4623e-09,
    -1.5845e-01,-2.0376e-02,-4.0081e-04, 2.2062e-04,-7.9179e-06,
    -1.6441e-07,-5.0004e-08, 8.0689e-10,-2.3813e-10,-2.4483e-10
  };
    
  /*static*/ constexpr double b_n[] = {
     0.0000e+00, 0.0000e+00,-1.1930e-02, 0.0000e+00, 9.8349e-03,
    -1.6861e-03, 0.0000e+00, 4.3338e-03, 6.1707e-03, 7.4635e-04,
     0.0000e+00, 3.5124e-03, 2.1967e-03, 4.2029e-04, 2.4476e-06,
     0.0000e+00, 4.1373e-04,-2.3281e-03, 2.7382e-04,-8.5220e-05,
     1.4204e-05, 0.0000e+00,-8.0076e-03, 4.5587e-05,-5.8053e-05,
    -1.1021e-05, 7.2338e-07,-1.9827e-07, 0.0000e+00,-3.9229e-03,
    -4.0697e-04,-1.6992e-04, 5.4705e-06,-4.4594e-06, 2.0121e-07,
    -7.7840e-08, 0.0000e+00,-3.2916e-03,-1.2302e-03,-6.5735e-06,
    -3.1840e-06,-8.9836e-07, 1.1870e-07,-5.8781e-09,-2.9124e-09,
     0.0000e+00, 1.0759e-02,-6.6074e-05,-4.0635e-05, 8.7141e-06,
     6.4567e-07,-4.4684e-08,-5.0293e-11, 2.7723e-10, 1.6903e-10
  };
    
  /*static*/ constexpr double a_e[] = { 
    -2.4104e-03, 1.1408e-04,-3.4621e-04, 1.6565e-03,-4.0620e-03,
    -6.8424e-03,-3.3718e-04, 7.3857e-03,-1.3324e-03,-1.5645e-03,
     4.6444e-03, 1.0296e-03, 3.6253e-03, 4.0329e-04, 3.1943e-04,
    -7.1992e-04, 4.8706e-03, 9.4300e-04, 2.0765e-04,-5.0987e-06,
    -7.1741e-06,-1.3131e-02, 2.9099e-04,-2.2509e-04, 2.6716e-04,
    -8.1815e-05, 8.4297e-06,-9.2378e-07,-5.8095e-04, 2.7501e-03,
     4.3659e-04,-8.2990e-06,-1.4808e-05, 2.2033e-06,-3.3215e-07,
     2.8858e-08, 9.9968e-03, 4.9291e-04, 3.3739e-05, 2.4696e-06,
    -8.1749e-06,-9.0052e-07, 2.0153e-07,-1.0271e-08, 1.8249e-09,
     3.0578e-03, 1.1229e-03,-1.9977e-04, 4.4581e-06,-7.6921e-06,
    -2.8308e-07, 1.0305e-07,-6.9026e-09, 1.5523e-10,-1.0395e-10
  };
    
  /*static*/ constexpr double b_e[] = {
     0.0000e+00, 0.0000e+00,-2.5396e-03, 0.0000e+00, 9.2146e-03,
    -7.5836e-03, 0.0000e+00, 1.2765e-02,-1.1436e-03, 1.7909e-04,
     0.0000e+00, 2.9318e-03,-6.8541e-04, 9.5775e-04, 2.4596e-05,
     0.0000e+00, 3.5662e-03,-1.3949e-03,-3.4597e-04,-5.8236e-05,
     5.6956e-06, 0.0000e+00,-5.0164e-04,-6.5585e-04, 1.1134e-05,
     2.3315e-05,-4.0521e-06,-4.1747e-07, 0.0000e+00, 5.1650e-04,
    -1.0483e-03, 5.8109e-06, 1.6406e-05,-1.6261e-06, 6.2992e-07,
     1.3134e-08, 0.0000e+00,-6.1449e-03,-3.2511e-04, 1.7646e-04,
     7.5326e-06,-1.1946e-06, 5.1217e-08, 2.4618e-08, 3.6290e-09,
     0.0000e+00, 3.6769e-03,-9.7683e-04,-3.2096e-07, 1.3860e-06,
    -6.2832e-09, 2.6918e-09, 2.5705e-09,-2.4401e-09,-3.7917e-11
  };
    
    
  // unit vector
  double x = cos(dlat) * cos(dlon);
  double y = cos(dlat) * sin(dlon);
  double z = sin(dlat);
    
  // Legendre polynomials
  double v[10][10],
         w[10][10];
  v[0][0] = 1e0;
  w[0][0] = 0e0;
  v[1][0] = z * v[0][0];
  w[1][0] = 0e0;

  int N, M;
  for (int n=1; n<nmax; n++) {
    N = n+1;
    v[n+1][0] = ((2*N-1)*z*v[n][0]-(N-1)*v[n-1][0]) / static_cast<double>(N);
    w[n+1][0] = 0e0;
  }
    
  for (int m=0; m<mmax; m++) {
    M = m+1;
    v[m+1][m+1] = static_cast<double>(2*M-1)*(x*v[m][m]-y*w[m][m]);
    w[m+1][m+1] = static_cast<double>(2*M-1)*(x*w[m][m]+y*v[m][m]);
    if (m<mmax-1) {
      v[m+2][m+1] = (2*M+1)*z*v[m+1][m+1];
      w[m+2][m+1] = (2*M+1)*z*w[m+1][m+1];
    }
    N = M+2;
    for (int n=m+2;n<nmax;n++) {
      v[n+1][m+1] = ((2*N-1)*z*v[n][m+1]-(N+M-1)*v[n-1][m+1]) 
        / static_cast<double>(N-M);
      w[n+1][m+1] = ((2*N-1)*z*w[n][m+1]-(N+M-1)*w[n-1][m+1]) 
        / static_cast<double>(N-M);
      N++;
    }
  }

  // Surface pressure on the geoid
  grn = 0e0;
  gre = 0e0;
  int i = 0;
  for (int n=0; n<=nmax; n++) {
    for (int m=0; m<=n; m++) {
      grn += (a_n[i]*v[n][m] + b_n[i]*w[n][m]);
      gre += (a_e[i]*v[n][m] + b_e[i]*w[n][m]);
      i++;
    }
  }
    
  // calculation of the asymmetric delay in m (Chen and Herring 1997)
  d = 1e0 / (sin(el)*tan(el)+0.0031e0)*(grn*cos(az)+gre*sin(az)); // mm
  d = d / 1000e0; // m 

  // Finished
  return 0;
}
