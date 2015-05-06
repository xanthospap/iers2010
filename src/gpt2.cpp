#include "iers2010.hpp"
#include <fstream>

#ifdef USE_EXTERNAL_CONSTS
    #include "gencon.hpp"
#endif

/** @brief  Define the path to the gpt2_5.grd file (including the filename)
 */
#define PATH_TO_GRD25_GRD "/usr/local/share/libiers10/gpt2_5.grd"

/** @brief         sign function or signum function; extracts the sign of a real number
 *  @param[in] val Template parameter
 *  @return        -1 if val < 0, <br>
 *                  0 if val = 0, <br>
 *                  1 if val > 0
 */
template <typename T> 
inline int sgn (T val) {
  return (T(0) < val) - (val < T(0));
}

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
 *                               -2 | Error reading gridfile
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
 * TODO Test some extreme cases
 */
int iers2010::gpt2 (const double& dmjd,const double* dlat,const double* dlon,
        const double* hell,const int& nstat,double* p,double* t,double* dt,
        double* e,double* ah,double* aw,double* undu,int it,const char* ifile)
{
    #ifdef USE_EXTERNAL_CONSTS
        constexpr double TWOPI   (D2PI);
        constexpr double PI      (DPI);
    #else
        //constexpr double TWOPI   (6.283185307179586476925287e0);
        constexpr double PI      (3.1415926535897932384626433e0);
        constexpr double TWOPI   ( 2.0e0 * PI );
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
    // lines in grid file
    constexpr int maxl = 2592;

    // Define factors for amplitudes
    double cosfy (0e0),coshy(0e0),sinfy(0e0),sinhy(0e0);
    if (it==1) { // constant parameters
        ;        // already initialized to zero
    } else {
        double fpy  (2e0 * TWOPI);
        double nom  (dmjd1 / 365.25e0);
        cosfy = cos (nom * TWOPI);
        coshy = cos (nom * fpy);
        sinfy = sin (nom * TWOPI);
        sinhy = sin (nom * fpy);
    }

    // Declare matrices to hold the grid
    double pgrid[maxl][5], tgrid[maxl][5], qgrid[maxl][5], dtgrid[maxl][5],
           ahgrid[maxl][5], awgrid[maxl][5], u[maxl], hs[maxl];

    // Read the external gridfile
    // The grid file was obtained from the website
    // http://acc.igs.org/tropo/gpt2_5.grd on 11/6/2012
    std::ifstream fin;
    if (!ifile)
        fin.open (PATH_TO_GRD25_GRD,std::ifstream::in);
    else
        fin.open (ifile,std::ifstream::in);
    if (!fin.is_open ())
        return -1;

    std::string s;
    // read the first comment line
    std::getline (fin,s);

    // Loop over grid points
    // WARNING Qgrid,  dTgrid, ahgrid and awgrid must be / 1000
    float dummy1,dummy2;
    for (int n=0;n<maxl;n++) {
        fin >> dummy1 >> dummy2
            >> pgrid[n][0] >> pgrid[n][1] >> pgrid[n][2] >> pgrid[n][3] >> pgrid[n][4]
            >> tgrid[n][0] >> tgrid[n][1] >> tgrid[n][2] >> tgrid[n][3] >> tgrid[n][4]
            >> qgrid[n][0] >> qgrid[n][1] >> qgrid[n][2] >> qgrid[n][3] >> qgrid[n][4]
            >>dtgrid[n][0] >>dtgrid[n][1] >>dtgrid[n][2] >>dtgrid[n][3] >>dtgrid[n][4]
            >> u[n] >> hs[n]
            >>ahgrid[n][0] >>ahgrid[n][1] >>ahgrid[n][2] >>ahgrid[n][3] >>ahgrid[n][4]
            >>awgrid[n][0] >>awgrid[n][1] >>awgrid[n][2] >>awgrid[n][3] >>awgrid[n][4]
            ;
        if (fin.fail ()) {
            fin.close ();
            return -2;
        }
    }

    // some variables to be used
    int indx[4];
    double undul[4],ql[4],dtl[4],tl[4],pl[4],ahl[4],awl[4];

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
        double diffpod ( (ppod - ( ipod*5e0 - 2.5e0))/5e0 );
        double difflon ( (plon - ( ilon*5e0 - 2.5e0))/5e0 );

        if (ipod==36)
            ipod -= 1;

        // get the number of the corresponding line
        indx[0] = (ipod-1)*72 + ilon;
        if (indx[0]>0) indx[0]--;

        // near the poles: nearest neighbour interpolation, otherwise: bilinear
        bool ibilinear (false);
        if ( (ppod>2.5e0) && (ppod<177.5e0) )
            ibilinear = true;

        // case of nearest neighbour
        if (!ibilinear) {

            int ix ( indx[0] );

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
            q /= 1000.0e0;

            // lapse rate of the temperature
            dt[k] = dtgrid[ix][0] +
                dtgrid[ix][1] * cosfy + dtgrid[ix][2] * sinfy +
                dtgrid[ix][3] * coshy + dtgrid[ix][4] * sinhy;
            dt[k] /= 1000.0e0;

            // station height - grid height
            double redh ( hgt - hs[ix] );

            // temperature at station height in Celsius
            t[k] = t0 + dt[k]*redh - 273.15e0;

            // temperature lapse rate in degrees / km
            dt[k] *= 1000e0;

            // virtual temperature in Kelvin
            double tv ( t0 * (1e0 + 0.6077e0*q) );

            double c ( GM * DMTR / (RG * tv) );

            // pressure in hPa
            p[k] = (p0*exp(-c*redh))/100e0;

            // water vapour pressure in hPa
            e[k] = (q*p[k])/(0.622e0 + 0.378e0*q);

            // hydrostatic coefficient ah 
            ah[k] = ahgrid[ix][0] +
                ahgrid[ix][1] * cosfy + ahgrid[ix][2] * sinfy +
                ahgrid[ix][3] * coshy + ahgrid[ix][4] * sinhy;
            ah[k] /= 1000.0e0;

            // wet coefficient aw
            aw[k] = awgrid[ix][0] +
                awgrid[ix][1] * cosfy + awgrid[ix][2] * sinfy +
                awgrid[ix][3] * coshy + awgrid[ix][4] * sinhy;
            aw[k] /= 1000.0e0;

        // bilinear interpolation
        } else {

            int ipod1 ( ipod + sgn<double>(diffpod) );
            int ilon1 ( ilon + sgn<double>(difflon) );
            if (ilon1==73)
                ilon1 = 1;
            if (!ilon)
                ilon1 = 72;

            // get the number of the line
            indx[1] = (ipod1 - 1)*72 + ilon;  // along same logtitude
            indx[2] = (ipod  - 1)*72 + ilon1; // along same polar distance
            indx[3] = (ipod1 - 1)*72 + ilon1; // diagonal
            if (indx[1]>0) indx[1]--;
            if (indx[2]>0) indx[2]--;
            if (indx[3]>0) indx[3]--;

            for (int l=0;l<4;l++) {

                int indxl = indx[l];

                // transforming ellipsoidial height to orthometric height:
                // Hortho = -N + Hell
                undul[l] = u[indxl];
                double hgt (hell[k]-undul[l]);

                // pressure, temperature at the heigtht of the grid
                double t0 = tgrid[indxl][0] +
                    tgrid[indxl][1] * cosfy + tgrid[indxl][2] * sinfy +
                    tgrid[indxl][3] * coshy + tgrid[indxl][4] * sinhy;
                double p0 = pgrid[indxl][0] +
                    pgrid[indxl][1] * cosfy + pgrid[indxl][2] * sinfy +
                    pgrid[indxl][3] * coshy + pgrid[indxl][4] * sinhy;

                // humidity
                ql[l] = qgrid[indxl][0] +
                    qgrid[indxl][1] * cosfy + qgrid[indxl][2] * sinfy +
                    qgrid[indxl][3] * coshy + qgrid[indxl][4] * sinhy;
                ql[l] /= 1000.0e0;

                // reduction = stationheight - gridheight
                double hs1  ( hs[indxl] );
                double redh ( hgt - hs1 );

                // lapse rate of the temperature in degree / m
                dtl[l] = dtgrid[indxl][0] +
                    dtgrid[indxl][1] * cosfy + dtgrid[indxl][2] * sinfy +
                    dtgrid[indxl][3] * coshy + dtgrid[indxl][4] * sinhy;
                dtl[l] /= 1000.0e0;

                // temperature reduction to station height
                tl[l] = t0 + dtl[l]*redh - 273.15e0;

                // virtual temperature
                double tv ( t0 * (1e0 + 0.6077e0 * ql[l]) );
                double c  ( GM*DMTR/(RG*tv) );

                // pressure in hPa
                pl[l] = (p0*exp(-c*redh))/100e0;

                // hydrostatic coefficient ah
                ahl[l] = ahgrid[indxl][0] +
                    ahgrid[indxl][1] * cosfy + ahgrid[indxl][2] * sinfy +
                    ahgrid[indxl][3] * coshy + ahgrid[indxl][4] * sinhy;
                ahl[l] /= 1000.0e0;

                // wet coefficient aw
                awl[l] = awgrid[indxl][0] +
                    awgrid[indxl][1] * cosfy + awgrid[indxl][2] * sinfy +
                    awgrid[indxl][3] * coshy + awgrid[indxl][4] * sinhy;
                awl[l] /= 1000.0e0;

            }

            double dnpod1 = std::abs (diffpod);  // distance nearer point
            double dnpod2 = 1.e0 - dnpod1;       // distance to distant point
            double dnlon1 = std::abs (difflon);
            double dnlon2 = 1.e0 - dnlon1;

            // pressure
            double r1 = dnpod2 * pl[0] + dnpod1 * pl[1];
            double r2 = dnpod2 * pl[2] + dnpod1 * pl[3];
            p[k] = dnlon2 * r1 + dnlon1 * r2;

            // temperature
            r1 = dnpod2 * tl[0] + dnpod1 * tl[1];
            r2 = dnpod2 * tl[2] + dnpod1 * tl[3];
            t[k] = dnlon2 * r1 + dnlon1 * r2;

            // temperature in degree per km
            r1 = dnpod2 * dtl[0] + dnpod1 * dtl[1];
            r2 = dnpod2 * dtl[2] + dnpod1 * dtl[3];
            dt[k] = (dnlon2 * r1 + dnlon1 * r2) * 1000.e0;

            // humidity
            r1 = dnpod2 * ql[0] + dnpod1 * ql[1];
            r2 = dnpod2 * ql[2] + dnpod1 * ql[3];
            double q ( dnlon2 * r1 + dnlon1 * r2 );
            e[k] = (q * p[k]) / (0.622e0 + 0.378e0 * q);

            // hydrostatic
            r1 = dnpod2 * ahl[0] + dnpod1 * ahl[1];
            r2 = dnpod2 * ahl[2] + dnpod1 * ahl[3];
            ah[k] = dnlon2 * r1 + dnlon1 * r2;

            // wet
            r1 = dnpod2 * awl[0] + dnpod1 * awl[1];
            r2 = dnpod2 * awl[2] + dnpod1 * awl[3];
            aw[k] = dnlon2 * r1 + dnlon1 * r2;

            // undulation
            r1 = dnpod2 * undul[0] + dnpod1 * undul[1];
            r2 = dnpod2 * undul[2] + dnpod1 * undul[3];
            undu[k] = dnlon2 * r1 + dnlon1 * r2;
        }
    }

    // Finished
    return 0;
}
