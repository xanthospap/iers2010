#include "iers2010.hpp"
#ifdef USE_EXTERNAL_CONSTS
  #include "gencon.hpp"
#endif

/**
 * @details This function evaluates the effects of zonal Earth tides on the
 *          rotation of the Earth.  The model used is a combination of Yoder
 *          et al. (1981) elastic body tide, Wahr and Bergen (1986) inelastic
 *          body tide, and Kantha et al. (1998) ocean tide models 
 *          as recommended by the IERS Conventions (2010).  Refer to
 *          Chapter 8 pp. 105 - 106.  The latest version of the model is located
 *          at http://tai.bipm.org/iers/convupdt/convupdt_c8.html.
 *          This function is a translation/wrapper for the fortran RG_ZONT2
 *          subroutine, found here : 
 *          http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  t      Julian centuries since J2000 (Note 1)
 * @param[out] dut    Effect on UT1 (Note 2)
 * @param[out] dlod   Effect on excess length of day (LOD) (Note 3)
 * @param[out] domega Effect on rotational speed (Note 4)
 * @return            An integer value always 0.
 * 
 * @note
 *  -# Though T is strictly TDB, it is usually more convenient to use
 *     TT, which makes no significant difference.  Julian centuries since
 *     J2000 is (JD - 2451545.0)/36525.
 *  -# The expression used is as adopted in IERS Conventions (2010).
 *     DUT is expressed in seconds and is double precision.
 *  -# The expression used is as adopted in IERS Conventions (2010).
 *     DLOD is the excess in LOD and is expressed in seconds per day
 *     and is double precision.  The phrase 'per day' is generally
 *     understood, so it has been omitted commonly in speech and
 *     literature.  
 *     See: Stephenson, F. R., Morrison, L. V., Whitrow, G. J., 1984,
 *     "Long-Term Changes in the Rotation of the Earth: 700 B. C. to
 *     A. D. 1980 [and Discussion]", Phil. Trans. Roy. Soc. of London.
 *     Series A, 313, pp. 47 - 70.
 *  -# The expression used is as adopted in IERS Conventions (2010).
 *     Rotational speed is expressed in radians per second and is
 *     double precision.
 *  -# Status:  Class 3 model
 * 
 * @version 20.12.2011
 * 
 * @cite iers2010, 
 *    Yoder, C. F., Williams, J. G., and Parke, M. E., (1981),
 *     "Tidal Variations of Earth Rotation," J. Geophys. Res., 86,
 *     pp. 881 - 891.
 *
 *     Wahr, J. and Bergen, Z., (1986), "The effects of mantle 
 *     anelasticity on nutations, Earth tides, and tidal variations
 *     in rotation rate," Geophys. J. Roy. astr. Soc., 87, pp. 633 - 668.
 *
 *     Kantha, L. H., Stewart, J. S., and Desai, S. D., (1998), "Long-
 *     period lunar fortnightly and monthly ocean tides," J. Geophys.
 *     Res., 103, pp. 12639 - 12647.
 *
 *     Gross, R. S., (2009), "Ocean tidal effects on Earth rotation,"
 *     J. Geodyn., 48(3-5), pp. 219 - 225.
 * 
 *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
 *     IERS Technical Note No. 36, BKG (2010)
 * 
 */
int
iers2010::rg_zont2(double t, double& dut, double& dlod, double& domega)
{

  // Set constants
#ifdef USE_EXTERNAL_CONSTS
  constexpr double TWOPI   (D2PI);
#else
  constexpr double TWOPI   (6.283185307179586476925287e0);
#endif

  /*  ----------------------
   *  Zonal Earth tide model
   *  ----------------------*/

  // Number of terms in the zonal Earth tide model  
  constexpr int nzont = 62;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   *  --------------------------------------------------
   *  Tables of multiples of arguments and coefficients
   *  --------------------------------------------------
   *  Luni-Solar argument multipliers
   *      l   l'  F   D OMEGA
   */
  constexpr int nfund[nzont][5] = {
    //  DATA ( ( NFUND(I,J), I=1,5 ), J= 1, 20 ) /
    { 1,  0,  2,  2,  2 },
    { 2,  0,  2,  0,  1 }, 
    { 2,  0,  2,  0,  2 }, 
    { 0,  0,  2,  2,  1 }, 
    { 0,  0,  2,  2,  2 }, 
    { 1,  0,  2,  0,  0 }, 
    { 1,  0,  2,  0,  1 }, 
    { 1,  0,  2,  0,  2 }, 
    { 3,  0,  0,  0,  0 }, 
    {-1,  0,  2,  2,  1 }, 
    {-1,  0,  2,  2,  2 }, 
    { 1,  0,  0,  2,  0 }, 
    { 2,  0,  2, -2,  2 }, 
    { 0,  1,  2,  0,  2 }, 
    { 0,  0,  2,  0,  0 }, 
    { 0,  0,  2,  0,  1 }, 
    { 0,  0,  2,  0,  2 }, 
    { 2,  0,  0,  0, -1 }, 
    { 2,  0,  0,  0,  0 }, 
    { 2,  0,  0,  0,  1 },             
    // DATA ( ( NFUND(I,J), I=1,5 ), J= 21, 40 ) /
    { 0, -1,  2,  0,  2 }, 
    { 0,  0,  0,  2, -1 }, 
    { 0,  0,  0,  2,  0 }, 
    { 0,  0,  0,  2,  1 }, 
    { 0, -1,  0,  2,  0 }, 
    { 1,  0,  2, -2,  1 }, 
    { 1,  0,  2, -2,  2 }, 
    { 1,  1,  0,  0,  0 }, 
    {-1,  0,  2,  0,  0 }, 
    {-1,  0,  2,  0,  1 }, 
    {-1,  0,  2,  0,  2 }, 
    { 1,  0,  0,  0, -1 }, 
    { 1,  0,  0,  0,  0 }, 
    { 1,  0,  0,  0,  1 }, 
    { 0,  0,  0,  1,  0 }, 
    { 1, -1,  0,  0,  0 }, 
    {-1,  0,  0,  2, -1 }, 
    {-1,  0,  0,  2,  0 }, 
    {-1,  0,  0,  2,  1 }, 
    { 1,  0, -2,  2, -1 },
    // DATA ( ( NFUND(I,J), I=1,5 ), J= 41, 60 ) /
    {-1, -1,  0,  2,  0 }, 
    { 0,  2,  2, -2,  2 }, 
    { 0,  1,  2, -2,  1 }, 
    { 0,  1,  2, -2,  2 }, 
    { 0,  0,  2, -2,  0 }, 
    { 0,  0,  2, -2,  1 }, 
    { 0,  0,  2, -2,  2 }, 
    { 0,  2,  0,  0,  0 }, 
    { 2,  0,  0, -2, -1 }, 
    { 2,  0,  0, -2,  0 }, 
    { 2,  0,  0, -2,  1 }, 
    { 0, -1,  2, -2,  1 }, 
    { 0,  1,  0,  0, -1 }, 
    { 0, -1,  2, -2,  2 }, 
    { 0,  1,  0,  0,  0 }, 
    { 0,  1,  0,  0,  1 }, 
    { 1,  0,  0, -1,  0 }, 
    { 2,  0, -2,  0,  0 }, 
    {-2,  0,  2,  0,  1 }, 
    {-1,  1,  0,  1,  0 },
    // DATA ( ( NFUND(I,J), I=1,5 ), J= 61, 62 ) /
    { 0,  0,  0,  0,  2 }, 
    { 0,  0,  0,  0,  1 }
  };

  /*
   *     Multiple of     DUT          DLOD              DOMEGA
   *            sin     cos      cos      sin       cos      sin
   */
  constexpr double tide[nzont][6] = {
    // DATA ( ( TIDE(I,J), I=1,6 ), J = 1,20 ) /
    {-0.0235e0,  0.0000e0, 0.2617e0, 0.0000e0, -0.2209e0, 0.0000e0},
    {-0.0404e0,  0.0000e0, 0.3706e0, 0.0000e0, -0.3128e0, 0.0000e0},
    {-0.0987e0,  0.0000e0, 0.9041e0, 0.0000e0, -0.7630e0, 0.0000e0},
    {-0.0508e0,  0.0000e0, 0.4499e0, 0.0000e0, -0.3797e0, 0.0000e0},
    {-0.1231e0,  0.0000e0, 1.0904e0, 0.0000e0, -0.9203e0, 0.0000e0},
    {-0.0385e0,  0.0000e0, 0.2659e0, 0.0000e0, -0.2244e0, 0.0000e0},
    {-0.4108e0,  0.0000e0, 2.8298e0, 0.0000e0, -2.3884e0, 0.0000e0},
    {-0.9926e0,  0.0000e0, 6.8291e0, 0.0000e0, -5.7637e0, 0.0000e0},
    {-0.0179e0,  0.0000e0, 0.1222e0, 0.0000e0, -0.1031e0, 0.0000e0},
    {-0.0818e0,  0.0000e0, 0.5384e0, 0.0000e0, -0.4544e0, 0.0000e0},
    {-0.1974e0,  0.0000e0, 1.2978e0, 0.0000e0, -1.0953e0, 0.0000e0},
    {-0.0761e0,  0.0000e0, 0.4976e0, 0.0000e0, -0.4200e0, 0.0000e0},
    { 0.0216e0,  0.0000e0,-0.1060e0, 0.0000e0,  0.0895e0, 0.0000e0},
    { 0.0254e0,  0.0000e0,-0.1211e0, 0.0000e0,  0.1022e0, 0.0000e0},
    {-0.2989e0,  0.0000e0, 1.3804e0, 0.0000e0, -1.1650e0, 0.0000e0},
    {-3.1873e0,  0.2010e0,14.6890e0, 0.9266e0,-12.3974e0,-0.7820e0},
    {-7.8468e0,  0.5320e0,36.0910e0, 2.4469e0,-30.4606e0,-2.0652e0},
    { 0.0216e0,  0.0000e0,-0.0988e0, 0.0000e0,  0.0834e0, 0.0000e0},
    {-0.3384e0,  0.0000e0, 1.5433e0, 0.0000e0, -1.3025e0, 0.0000e0},
    { 0.0179e0,  0.0000e0,-0.0813e0, 0.0000e0,  0.0686e0, 0.0000e0},
    // DATA ( ( TIDE(I,J), I=1,6 ), J = 21,40 ) /
    {-0.0244e0,  0.0000e0, 0.1082e0, 0.0000e0, -0.0913e0, 0.0000e0},
    { 0.0470e0,  0.0000e0,-0.2004e0, 0.0000e0,  0.1692e0, 0.0000e0},
    {-0.7341e0,  0.0000e0, 3.1240e0, 0.0000e0, -2.6367e0, 0.0000e0},
    {-0.0526e0,  0.0000e0, 0.2235e0, 0.0000e0, -0.1886e0, 0.0000e0},
    {-0.0508e0,  0.0000e0, 0.2073e0, 0.0000e0, -0.1749e0, 0.0000e0},
    { 0.0498e0,  0.0000e0,-0.1312e0, 0.0000e0,  0.1107e0, 0.0000e0},
    { 0.1006e0,  0.0000e0,-0.2640e0, 0.0000e0,  0.2228e0, 0.0000e0},
    { 0.0395e0,  0.0000e0,-0.0968e0, 0.0000e0,  0.0817e0, 0.0000e0},
    { 0.0470e0,  0.0000e0,-0.1099e0, 0.0000e0,  0.0927e0, 0.0000e0},
    { 0.1767e0,  0.0000e0,-0.4115e0, 0.0000e0,  0.3473e0, 0.0000e0},
    { 0.4352e0,  0.0000e0,-1.0093e0, 0.0000e0,  0.8519e0, 0.0000e0},
    { 0.5339e0,  0.0000e0,-1.2224e0, 0.0000e0,  1.0317e0, 0.0000e0},
    {-8.4046e0,  0.2500e0,19.1647e0, 0.5701e0,-16.1749e0,-0.4811e0},
    { 0.5443e0,  0.0000e0,-1.2360e0, 0.0000e0,  1.0432e0, 0.0000e0},
    { 0.0470e0,  0.0000e0,-0.1000e0, 0.0000e0,  0.0844e0, 0.0000e0},
    {-0.0555e0,  0.0000e0, 0.1169e0, 0.0000e0, -0.0987e0, 0.0000e0},
    { 0.1175e0,  0.0000e0,-0.2332e0, 0.0000e0,  0.1968e0, 0.0000e0},
    {-1.8236e0,  0.0000e0, 3.6018e0, 0.0000e0, -3.0399e0, 0.0000e0},
    { 0.1316e0,  0.0000e0,-0.2587e0, 0.0000e0,  0.2183e0, 0.0000e0},
    { 0.0179e0,  0.0000e0,-0.0344e0, 0.0000e0,  0.0290e0, 0.0000e0},
    // DATA ( ( TIDE(I,J), I=1,6 ), J = 41,60 ) /
    {-0.0855e0,  0.0000e0, 0.1542e0, 0.0000e0, -0.1302e0, 0.0000e0},
    {-0.0573e0,  0.0000e0, 0.0395e0, 0.0000e0, -0.0333e0, 0.0000e0},
    { 0.0329e0,  0.0000e0,-0.0173e0, 0.0000e0,  0.0146e0, 0.0000e0},
    {-1.8847e0,  0.0000e0, 0.9726e0, 0.0000e0, -0.8209e0, 0.0000e0},
    { 0.2510e0,  0.0000e0,-0.0910e0, 0.0000e0,  0.0768e0, 0.0000e0},
    { 1.1703e0,  0.0000e0,-0.4135e0, 0.0000e0,  0.3490e0, 0.0000e0},
    {-49.7174e0, 0.4330e0,17.1056e0, 0.1490e0,-14.4370e0,-0.1257e0},
    {-0.1936e0,  0.0000e0, 0.0666e0, 0.0000e0, -0.0562e0, 0.0000e0},
    { 0.0489e0,  0.0000e0,-0.0154e0, 0.0000e0,  0.0130e0, 0.0000e0},
    {-0.5471e0,  0.0000e0, 0.1670e0, 0.0000e0, -0.1409e0, 0.0000e0},
    { 0.0367e0,  0.0000e0,-0.0108e0, 0.0000e0,  0.0092e0, 0.0000e0},
    {-0.0451e0,  0.0000e0, 0.0082e0, 0.0000e0, -0.0069e0, 0.0000e0},
    { 0.0921e0,  0.0000e0,-0.0167e0, 0.0000e0,  0.0141e0, 0.0000e0},
    { 0.8281e0,  0.0000e0,-0.1425e0, 0.0000e0,  0.1202e0, 0.0000e0},
    { -15.8887e0,0.1530e0, 2.7332e0, 0.0263e0, -2.3068e0,-0.0222e0},
    {-0.1382e0,  0.0000e0, 0.0225e0, 0.0000e0, -0.0190e0, 0.0000e0},
    { 0.0348e0,  0.0000e0,-0.0053e0, 0.0000e0,  0.0045e0, 0.0000e0},
    {-0.1372e0,  0.0000e0,-0.0079e0, 0.0000e0,  0.0066e0, 0.0000e0},
    { 0.4211e0,  0.0000e0,-0.0203e0, 0.0000e0,  0.0171e0, 0.0000e0},
    {-0.0404e0,  0.0000e0, 0.0008e0, 0.0000e0, -0.0007e0, 0.0000e0},
    // DATA ( ( TIDE(I,J), I=1,6 ), J = 61,62 ) /
    { 7.8998e0,  0.0000e0, 0.1460e0, 0.0000e0, -0.1232e0, 0.0000e0},
    {-1617.2681e0,0.0000e0,-14.9471e0,0.0000e0, 12.6153e0, 0.0000e0}
  };

  // Computation of fundamental arguments
  double fargs[5]; /* (l, lp, f, d, om) */
  iers2010::fundarg(t, fargs);
  double l  = fargs[0];
  double lp = fargs[1];
  double f  = fargs[2];
  double d  = fargs[3];
  double om = fargs[4];

  // Set initial values to zero.
  dut    = 0e0;
  dlod   = 0e0;
  domega = 0e0;

  //  Sum zonal tide terms.
  double arg, dsinarg, dcosarg;
  for (int i=0; i<nzont; i++) {
    // Formation of multiples of arguments.
    arg = std::fmod(
        (double)( nfund[i][0] ) * l 
        + (double)( nfund[i][1] ) * lp 
        + (double)( nfund[i][2] ) * f 
        + (double)( nfund[i][3] ) * d 
        + (double)( nfund[i][4] ) * om, TWOPI );
    while (arg < 0e0) { arg += TWOPI; }

    // Evaluate zonal tidal terms.
    dsinarg = std::sin(arg);
    dcosarg = std::cos(arg);
    dut    += tide[i][0] *dsinarg + tide[i][1] *dcosarg;
    dlod   += tide[i][2] *dcosarg + tide[i][3] *dsinarg;
    domega += tide[i][4] *dcosarg + tide[i][5] *dsinarg;
  }

  // Rescale corrections so that they are in units involving seconds.
  dut    *= 1e-4;
  dlod   *= 1e-5;
  domega *= 1e-14;

  //  Finsihed
  return 0;
}
