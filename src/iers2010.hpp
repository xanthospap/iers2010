#include <cmath>

namespace iers2010
{
    /// Compute the lunisolar fundamental arguments.
    int fundarg(double jc, double* fargs);

    /// Compute the diurnal lunisolar effect on polar motion.
    int pmsdnut2(double mjd, double& dx, double& dy);

    /// Compute the subdiurnal librations in UT1.
    int utlibr(double mjd, double& dut1, double& dlod);

    /// Compute corrections to the coordinates of the CIP to account for 
    /// Free Core Nutation.
    int fcnnut(double mjd, double& x, double& y, double& dx, double& dy);

    /// Compute the angular argument which depends on time for 11 tidal
    /// argument calculations.
    int arg2(int year, double day, double* angular_arg);

    /// Functions needed by dehanttideinel.
    namespace dhtide {
        
        /// Gregorian Calendar to Julian Date.
        int cal2jd(int iy, int im, int id, double& djm0, double& djm);

        /// For a given UTC date, calculate delta(AT) = TAI-UTC.
        int dat(int iy, int im, int id, double fd, double& deltat);

        /// Compute out-of-phase corrections induced by mantle anelasticity in
        /// the diurnal band. 
        void st1idiu(const double* xsta,const double* xsun, const double* xmon,
            double fac2sun, double fac2mon, double* xcorsta);

        /// Compute out-of-phase corrections induced by mantle anelasticity in 
        /// the semi-diurnal band.
        void st1isem(const double* xsta, const double* xsun, const double* xmon,
            double fac2sun, double fac2mon, double* xcorsta);

        /// Compute the corrections induced by the latitude dependence given
        /// by L^1 in Mathews et al. 1991
        void st1l1(const double* xsta,const double* xsun, const double* xmon,
            double fac2sun, double fac2mon, double* xcorsta);

        /// Compute in-phase and out-of-phase corrections induced by mantle 
        /// anelasticity in the diurnal band.
        void step2diu(const double* xsta, double fhr, double t, double* xcorsta);

        /// Compute the in-phase and out-of-phase corrections induced by mantle
        /// anelasticity in the long period band.
        void step2lon(const double* xsta, double t, double* xcorsta);
    
    } // dhtide

    /// Compute tidal corrections of station displacements caused by lunar and
    /// solar gravitational attraction.
    int dehanttideinel(const double* xsta,const double* xsun, const double* xmon,
        int yr, int month, int day, double fhr, double* dxtide);

    namespace hisp {

        /// Compute the difference between Epheremis Time and (UTC).
        double etutc(double year);

        /// Compute the frequency and phase of a tidal constituent given its 
        /// Doodson number.
        int tdfrph(const int idood[6], const int itm[5], double& freq, double& phase);

        /// Find an array s for the spline interpolation.
        int spline(int n, const double* x, const double* u, double* s, double* a);

        /// Perform cubic spline interpolation of a given function sampled at
        /// unequally spaced intervals.
        double eval(double y, int nn, const double* x, const double* u, 
            const double* s);

        /// Fill in data x, with sines and cosines of frequencies om.
        int recurs(double* x, int n, const double* hc, int nf, const double* om,
            double* scr);

        /// Sort array x and store index keys in k.
        int shells(double* x, int* k, int n);

        int admint(const double* ampin, const int idtin[][6], const double* phin,
            double* amp, double* f, double* p, int nin, int& nout, const int itm[5]);
 
    } // hisp

    namespace oeop {

        /// Compute the time dependent part of second degree diurnal and 
        /// semidiurnal tidal potential.
        int cnmtx(double dmjd, double* h);

    } // oeop

    /// Compute the diurnal and semi-diurnal variations in Earth Orientation
    /// Parameters from ocean tides.
    int ortho_eop(double time, double& dx, double& dy, double& dut1);

    /// Evaluate the effects of zonal Earth tides on the rotation of the Earth.
    int rg_zont2(double t, double& dut, double& dlod, double& domega);

    /// Compute the global total FCULa mapping function.
    double fcul_a(double dlat, double dhgt, double t, double elev);

    /// Computes the global total FCULb mapping function.
    double fcul_b(double dlat, double dhgt, double doy, double elev);

    /// Determines the total zenith delay following (Mendes and Pavlis, 2004).
    int fculzd_hpa(double dlat, double dhgt, double pres, double wvp, 
        double lambda, double& f_ztd, double& f_zhd, double& f_zwd);

    /// Compute the Global Mapping Functions (GMF).
    int gmf(double mjd, double lat, double lon, double hgt, double zd,
        double& gmfh, double& gmfw);

    /// Compute the Vienna Mapping Functions 1 (VMF1), to be used with "a" 
    /// coefficients computed for a given site.
    int vmf1(double ah, double aw, double dmjd, double dlat, double zd,
        double& vmf1h, double& vmf1w);

    /// Compute the Vienna Mapping Functions 1 (VMF1), with height corrections,
    /// to be used with "a" coefficients computed for a grid.
    int vmf1_ht(double ah, double aw, double dmjd, double dlat, double ht,
        double zd, double& vmf1h, double& vmf1w);

    /// Compute the Global Pressure and Temperature (GPT), based on spherical 
    /// harmonics up to degree and order 9. 
    int gpt(double dmjd, double dlat, double dlon, double dhgt, double& pres,
        double& temp, double& undu);

} // iers2010
