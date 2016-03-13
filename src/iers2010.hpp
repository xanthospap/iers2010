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

} // iers2010
