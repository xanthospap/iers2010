#include <cmath>

namespace iers2010
{
    /// Compute the lunisolar fundamental arguments.
    int fundarg(double jc, double* fargs);

    /// Evaluate polar motion due to tidal gravitation.
    int pmsdnut2(double mjd, double& dx, double& dy);

    /// Evaluates subdiurnal libration (dut1, dlod) due to tidal gravitation.
    int utlibr(double mjd, double& dut1, double& dlod);

} // iers2010
