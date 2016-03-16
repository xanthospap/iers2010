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

} // iers2010
