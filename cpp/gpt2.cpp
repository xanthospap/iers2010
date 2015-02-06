#include "iers2010.hpp"

#ifdef USE_EXTERNAL_CONSTS
    #include "gencon.hpp"
#endif

/**
 * @details  This function determines pressure, temperature, temperature lapse
 *           rate, water vapour pressure, hydrostatic and wet mapping function 
 *           coefficients ah and aw, and geoid undulation for specific sites 
 *           near the Earth surface. It is based on a 5 x 5 degree external 
 *           grid file ('gpt2_5.grd') with mean values as well as sine and 
 *           cosine amplitudes for the annual and semiannual variation of the 
 *           coefficients.
 *           This function is a translation/wrapper for the fortran GPT2
 *           subroutine, found here : 
 *           http://maia.usno.navy.mil/conv2010/software.html
 * 
 * @param[in]  dmjd  Modified Julian Date
 * @param[in]  dlat  Ellipsoidal latitude given in radians [-pi/2:+pi/2] 
 *                   (vector)
 * @param[in]  dlon  Longitude given in radians [-pi:pi] or [0:2pi] (vector)
 * @param[in]  hell  Ellipsoidal height in meters (vector)
 * @param[in]  nstat Number of stations in DLAT, DLON, and HELL (i.e. size of
 *                   arrays)
 * @param[in]  it    An integer, deonting:<br>
 *                   case 1 : no time variation but static quantities<br>
 *                   case 0 : with time variation (annual and semiannual terms)
 * @param[in]  ifile The name of the input gridfile to read values from. By
 *                   default, it is set to 'gpt2_5.grd'.
 * @param[out] p     Pressure given in hPa (vector of length nstat)
 * @param[out] t     Temperature in degrees Celsius (vector of length nstat)
 * @param[out] dt    Temperature lapse rate in degrees per km (vector of length
 *                   nstat)
 * @param[out] e     Water vapour pressure in hPa (vector of length nstat)
 * @param[out] ah    Hydrostatic mapping function coefficient at zero height 
 *                   (VMF1) (vector of length NSTAT)
 * @param[out] aw    Wet mapping function coefficient (VMF1) (vector of length 
 *                   NSTAT)
 * @param[out] undu  Geoid undulation in meters (vector of length nstat)
 * @return           An integer which can be:
 *                   Returned Value | Status
 *                   ---------------|-----------------------------------
 *                               -1 | Unable to open input gridfile
 * 
 * @note
 *    -# The hydrostatic mapping function coefficients have to be used with the
 *       height dependent Vienna Mapping Function 1 (vmf_ht.f) because the
 *       coefficients refer to zero height.
 *    -# Status: Class 1 model
 *
 * @verbatim
 * Test cases:
 * Example 1 (Vienna, 2 August 2012, with time variation):
 * dmjd = 56141.d0
 * dlat(1) = 48.20d0*pi/180.d0
 * dlon(1) = 16.37d0*pi/180.d0
 * hell(1) = 156.d0
 * nstat = 1
 * it = 0
 * 
 * output:
 * p = 1002.56 hPa
 * T = 22.12 deg Celsius
 * dT = -6.53 deg / km
 * e = 15.63 hPa
 * ah = 0.0012647
 * aw = 0.0005726
 * undu = 44.06 m
 * 
 * Example 2 (Vienna, 2 August 2012, without time variation, 
 * i.e. constant values):
 * 
 * dmjd = 56141.d0
 * dlat(1) = 48.20d0*pi/180.d0
 * dlon(1) = 16.37d0*pi/180.d0
 * hell(1) = 156.d0
 * nstat = 1
 * it = 1
 *
 * output:
 * p = 1003.49 hPa
 * T = 11.95 deg Celsius
 * dT = -5.47 deg / km
 * e = 9.58 hPa
 * ah = 0.0012395
 * aw = 0.0005560
 * undu = 44.06 m
 * @endverbatim
 * 
 * @version 2013 May 31
 * 
 * @cite iers2010
 *     Lagler, K., Schindelegger, M., Boehm, J., Krasna, H., and Nilsson, T.,
 *     (2013), "GPT2: Empirical slant delay model for radio space geodetic
 *     techniques," Geophys. Res. Lett., Vol. 40, pp. 1069-1073, DOI:
 *     10.1002/grl.50288. 
 * 
 */
int iers2010::gpt2 (const double& dmjd,const double* dlat,const double* dlon,
        const double* hell,const int& nstat,double* p,double* t,double* dt,
        double* e,double* ah,double* aw,double* undu,int it=0,const char* ifile)
{
    #ifdef USE_EXTERNAL_CONSTS
        constexpr double TWOPI   (D2PI);
        constexpr double PI      (DPI);
    #else
        constexpr double TWOPI   (6.283185307179586476925287e0);
        constexpr double PI      (3.1415926535e0);
    #endif

    // quick return, if nstat less than 1
    if (nstat<1)
        return 0;

    // Define the mean gravity in m/s**2
    constexpr double GM = 9.80665e0;
    // Define the molar mass of dry air in kg/mol
    constexpr double DMTR = 28.965e-3;
    // Universal gas constant in J/K/mol
    constexpr double RG = 8.3143e0;
    // Change the reference epoch to January 1 2000
    double dmjd1 = dmjd - 51544.5e0;
    // (180/PI) is used all the time !
    constexpr double d2r (180e0/PI);

    // Define factors for amplitudes
    double cosfy (0e0),coshy(0e0),sinfy(0e0),sinhy(0e0);
    if (it==1) { // constant parameters
        ;        // already initialized to zero
    } else {
        double tpy (365.25e0 * TWOPI);
        double fpy (365.25e0 * 2e0 * TWOPI);
        cosfy = cos (dmjd1 / tpy);
        coshy = cos (dmjd1 / fpy);
        sinfy = sin (dmjd1 / tpy);
        sinhy = sin (dmjd1 / fpy);
    }

    // Read the external gridfile
    // The grid file was obtained from the website
    // http://acc.igs.org/tropo/gpt2_5.grd on 11/6/2012
    std::ifstream fin;
    if (ifile==nullptr)
        fin.open ("gpt2_5.grd",std::ifstream::in);
    else
        fin.open (ifile,std::ifstream::in);
    if (!fin.is_open ())
        return -1

    char line[256];
    // read the first comment line
    fin.getline (line,256);

    // Loop over grid points
    fin.getline (line,256);
    for (int n=0;n<2592;n++) {
        pgrid[n] = 
        tgrid[n]
        qgrid[n]
        dtgrid[n]
        u[n]
        hs[n]
        ahgrid[n]
        awgrid[n]
        fin.getline (line,256);
    }

    // Loop over stations
    for (int k=0;k<nstat;k++) {
        
        // only positive longitude in degrees
        double plon ( (dlon[k]<0e0) ? ((dlon[k]+TWOPI)/d2r) : (dlon[k]*d2r) );
        // transform to polar distance in degrees
        double ppod ( (-dlat[k] + PI/2e0) * d2r );

        // find the index (line in the grid file) of the nearest point
        int ipod = floor ( (ppod+5e0)/5e0 );
        int ilon = floor ( (plon+5e0)/5e0 );

        // normalized (to one) differences, can be positive or negative
        double diffpod ( (ppod - (ipod*5e0 - 2.5e0))/5e0 );
        double difflon ( (plon - (ilon*5e0 - 2.5e0))/5e0 );

        if (ipod==37)
            ipod -= 1;

        // get the number of the corresponding line
        indx[0] = (ipod-1)*72 + ilon;

        // near the poles: nearest neighbour interpolation, otherwise: bilinear
        bool ibilinear (false);
        if ( (ppod>2.5e0) && (ppod<177.5e0) )
            ibilinear = true;

        // case of nearest neighbour
        if (!ibilinear) {

            int ix (indx[0]);

            // transforming ellipsoidial height to orthometric height
            undu[k] = u[ix];
            double hgt ( hell[k] - undu[k] );

            // pressure, temperature at the height of the grid
            double t0 = tgrid[ix][0] +
                tgrid[ix][1] * cosfy + tgrid[ix][2] * sinfy +
                tgrid[ix][3] * coshy + tgrid[ix][4] * sinhy;
            double p0 = pgrid[ix][0] +
                pgrid[ix][1] * cosfy + pgrid[ix][2] * sinfy +
                pgrid[ix][3] * coshy + pgrid[ix][4] * sinhy;

            // specify humidity
            double q = qgrid[ix][0] +
                qgrid[ix][1] * cosfy + qgrid[ix][2] * sinfy +
                qgrid[ix][3] * coshy + qgrid[ix][4] * sinhy;

            // lapse rate of the temperature
            dt[k] = dtgrid[ix][0] +
                dtgrid[ix][1] * cosfy + dtgrid[ix][2] * sinfy +
                dtgrid[ix][3] * coshy + dtgrid[ix][4] * sinhy;

            // station height - grid height
            double redh ( hgt - hs[ix] );

            // temperature at station height in Celsius
            t[k] = t0 + dt[k]*redh - 273.15e0;

            // temperature lapse rate in degrees / km
            dt[k] *= 1000e0;

            // virtual temperature in Kelvin
            double tv ( t0 * (1e0 + 0.6077e0*q) );

            double c ( GM * DTMR / (RG * tv) );

            // pressure in hPa
            p[k] = (p0*exp(-c*redh))/100e0;

            // water vapour pressure in hPa
            e[k] = (q*p[k])/(0.622e0 + 0.378e0*q);

            // hydrostatic coefficient ah 
            ah[k] = ahgrid[ix][0] +
                ahgrid[ix][1] * cosfy + ahgrid[ix][2] * sinfy +
                ahgrid[ix][3] * coshy + ahgrid[ix][4] * sinhy;

            // wet coefficient aw
            aw[k] = awgrid[ix][0] +
                awgrid[ix][1] * cosfy + awgrid[ix][2] * sinfy +
                awgrid[ix][3] * coshy + awgrid[ix][4] * sinhy;

        // bilinear interpolation
        } else {

        }

    
    // Finished
    return 0;
}
