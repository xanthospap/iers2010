#include "iers2010.hpp"
#include "matvec.hpp"

using dso::Vector3;

// Set constants
constexpr const double DEG2RAD = iers2010::D2PI / 360e0;

/// @brief Structure to hold a number of angles used to compute the Step 2
///        tidal displacement for a point on earth. See function step2diu and
///        step2lon. Instead of computing them within these functions, they are
///        computed once, and stored in a Step2Angles instance.
struct Step2Angles {
  Step2Angles(double t, double fhr) noexcept {
    // Compute the phase angles in degrees.
    s = 218.31664563e0 +
        (481267.88194e0 + (-0.0014663889e0 + (0.00000185139e0) * t) * t) * t;

    tau = fhr * 15e0 + 280.4606184e0 +
          (36000.7700536e0 + (0.00038793e0 + (-0.0000000258e0) * t) * t) * t +
          (-s);

    const double pr =
        (1.396971278e0 +
         (0.000308889e0 + (0.000000021e0 + (0.000000007e0) * t) * t) * t) *
        t;

    s += pr;

    h = 280.46645e0 +
        (36000.7697489e0 +
         (0.00030322222e0 + (0.000000020e0 + (-0.00000000654e0) * t) * t) * t) *
            t;

    p = 83.35324312e0 +
        (4069.01363525e0 +
         (-0.01032172222e0 + (-0.0000124991e0 + (0.00000005263e0) * t) * t) *
             t) *
            t;

    zns = 234.95544499e0 +
          (1934.13626197e0 +
           (-0.00207561111e0 + (-0.00000213944e0 + (0.00000001650e0) * t) * t) *
               t) *
              t;

    ps = 282.93734098e0 +
         (1.71945766667e0 +
          (0.00045688889e0 + (-0.00000001778e0 + (-0.00000000334e0) * t) * t) *
              t) *
             t;

    // Reduce angles to between the range 0 and 360.
    s = std::fmod(s, 360e0);
    tau = std::fmod(tau, 360e0);
    h = std::fmod(h, 360e0);
    p = std::fmod(p, 360e0);
    zns = std::fmod(zns, 360e0);
    ps = std::fmod(ps, 360e0);
  }
  double s, tau, h, p, zns, ps;
}; // Step2Angles

/// @brief A list of auxiliary quantities used to compute Step-1 and Step-2
///        tidal displacement. Originaly, they are performed within each of the
///        step-1 functions (st1l1, st1isem, st1diu) and step-1 functions (
///        step2diu and step2lon) but we can only compute them once and store
///        the in a TideAux instance.
struct TideAux {
  TideAux(const Vector3 *xsta, const Vector3 *xsun, const Vector3 *xmon,
          double fsun, double fmon) noexcept {
    rsta = xsta->norm();
    sinphi = xsta->z() / rsta;
    cosphi = std::sqrt(xsta->x() * xsta->x() + xsta->y() * xsta->y()) / rsta;
    sinla = xsta->y() / cosphi / rsta;
    cosla = xsta->x() / cosphi / rsta;
    rsun = xsun->norm();
    rmon = xmon->norm();
    fac2sun = fsun;
    fac2mon = fmon;
  };

  double rsta;
  double rsun;
  double rmon;
  double sinphi;
  double cosphi;
  double sinla;
  double cosla;
  double fac2sun;
  double fac2mon;
}; // TideAux

/// @details This function gives the corrections induced by the latitude
///          dependence given by L^1 in Mathews et al. 1991 (See References).
///          This function is a translation/wrapper for the fortran ST1L1
///          subroutine, found here :
///          http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  xsta    Geocentric position of the IGS station (Note 1)
/// @param[in]  xsun    Geocentric position of the Sun (Note 2)
/// @param[in]  xmon    Geocentric position of the Moon (Note 2)
/// @param[in]  fac2sun Degree 2 TGP factor for the Sun (Note 3)
/// @param[in]  fac2mon Degree 2 TGP factor for the Moon (Note 3)
/// @return  xcorsta Out of phase station corrections for semi-diurnal band
///
/// @note
///    -# The IGS station is in ITRF co-rotating frame. All coordinates are
///       expressed in meters, as arrays, i.e. [x,y,z].
///    -# The position is in Earth Centered Earth Fixed (ECEF) frame.  All
///       coordinates are expressed in meters, as arrays, i.e. [x,y,z].
///    -# The expressions are computed in the main program. TGP is the tide
///       generated potential. The units are inverse meters.
///    -# Status: Class 1
///    -# This fucnction is part of the package dehanttideinel, see
///       ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/
///
/// @version 31.07.2009
///
/// @cite iers2010,
///       Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, "Tidal station
///       displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477,
///       Mathews, P. M., Buffett, B. A., Herring, T. A., Shapiro, I. I.,
///       1991b, Forced nutations of the Earth: Influence of inner core
///       Dynamics 2. Numerical results and comparisons, J. Geophys. Res.,
///       96, 8243-8257
Vector3 st1l1(/*const Vector3 &xsta,*/ const Vector3 &xsun, const Vector3 &xmon,
              const TideAux &aux) noexcept {
  constexpr const double l1d = 0.0012e0;
  constexpr const double l1sd = 0.0024e0;

  // Compute the normalized position vector of the station.
  // const double rsta = aux.rsta;
  const double sinphi = aux.sinphi;
  const double sinphi2 = sinphi * sinphi;
  const double cosphi = aux.cosphi;
  const double cosphi2 = cosphi * cosphi;
  const double sinla = aux.sinla;
  const double cosla = aux.cosla;

  // Compute the normalized position vector of the Moon.
  const double rmon2 = aux.rmon * aux.rmon;

  // Compute the normalized position vector of the Sun.
  const double rsun2 = aux.rsun * aux.rsun;

  // Compute the station corrections for the diurnal band.
  double l1 = l1d;
  const double fac2sun = aux.fac2sun;
  const double fac2mon = aux.fac2mon;
  double dnsun = -l1 * sinphi2 * fac2sun * xsun.z() *
                 (xsun.x() * cosla + xsun.y() * sinla) / rsun2;
  double dnmon = -l1 * sinphi2 * fac2mon * xmon.z() *
                 (xmon.x() * cosla + xmon.y() * sinla) / rmon2;
  double desun = l1 * sinphi * (cosphi2 - sinphi2) * fac2sun * xsun.z() *
                 (xsun.x() * sinla - xsun.y() * cosla) / rsun2;
  double demon = l1 * sinphi * (cosphi2 - sinphi2) * fac2mon * xmon.z() *
                 (xmon.x() * sinla - xmon.y() * cosla) / rmon2;

  double de = 3e0 * (desun + demon);
  double dn = 3e0 * (dnsun + dnmon);

  Vector3 xcorsta{{-de * sinla - dn * sinphi * cosla,
                   de * cosla - dn * sinphi * sinla, dn * cosphi}};

  // Compute the station corrections for the semi-diurnal band.
  l1 = l1sd;
  const double costwola = cosla * cosla - sinla * sinla;
  const double sintwola = 2e0 * cosla * sinla;
  const double xsun_x2 = std::pow(xsun.x(), 2);
  const double xsun_y2 = std::pow(xsun.y(), 2);
  const double xmon_x2 = std::pow(xmon.x(), 2);
  const double xmon_y2 = std::pow(xmon.y(), 2);

  dnsun =
      -l1 / 2e0 * sinphi * cosphi * fac2sun *
      ((xsun_x2 - xsun_y2) * costwola + 2e0 * xsun.x() * xsun.y() * sintwola) /
      rsun2;

  dnmon =
      -l1 / 2e0 * sinphi * cosphi * fac2mon *
      ((xmon_x2 - xmon_y2) * costwola + 2e0 * xmon.x() * xmon.y() * sintwola) /
      rmon2;

  desun =
      -l1 / 2e0 * sinphi2 * cosphi * fac2sun *
      ((xsun_x2 - xsun_y2) * sintwola - 2e0 * xsun.x() * xsun.y() * costwola) /
      rsun2;

  demon =
      -l1 / 2e0 * sinphi2 * cosphi * fac2mon *
      ((xmon_x2 - xmon_y2) * sintwola - 2e0 * xmon.x() * xmon.y() * costwola) /
      rmon2;

  de = 3e0 * (desun + demon);
  dn = 3e0 * (dnsun + dnmon);

  xcorsta.data[0] += (-de * sinla - dn * sinphi * cosla);
  xcorsta.data[1] += (de * cosla - dn * sinphi * sinla);
  xcorsta.data[2] += (dn * cosphi);

  // Finished
  return xcorsta;
}

/// @details This function gives the out-of-phase corrections induced by
///          mantle anelasticity in the semi-diurnal band.
///          This function is a translation/wrapper for the fortran ST1ISEM
///          subroutine, found here :
///          http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  xsta    Geocentric position of the IGS station (Note 1)
/// @param[in]  xsun    Geocentric position of the Sun (Note 2)
/// @param[in]  xmon    Geocentric position of the Moon (Note 2)
/// @param[in]  fac2sun Degree 2 TGP factor for the Sun (Note 3)
/// @param[in]  fac2mon Degree 2 TGP factor for the Moon (Note 3)
/// @return xcorsta Out of phase station corrections for semi-diurnal band
///
/// @note
///     -# The IGS station is in ITRF co-rotating frame. All coordinates are
///        expressed in meters, as arrays, i.e. [x,y,z].
///     -# The position is in Earth Centered Earth Fixed (ECEF) frame.  All
///        coordinates are expressed in meters, as arrays, i.e. [x,y,z].
///     -# The expressions are computed in the main program. TGP is the tide
///        generated potential. The units are inverse meters.
///     -# Status: Class 1
///     -# This fucnction is part of the package dehanttideinel, see
///        ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/
///
/// @version 31.07.2009
///
/// @cite iers2010,
///       Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, "Tidal station
///       displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
Vector3 st1isem(/*const Vector3 &xsta,*/ const Vector3 &xsun,
                const Vector3 &xmon, const TideAux &aux) noexcept {

  constexpr const double dhi = -0.0022e0;
  constexpr const double dli = -0.0007e0;

  // Compute the normalized position vector of the IGS station.
  // const double rsta = aux.rsta;

  const double sinphi = aux.sinphi;
  const double cosphi = aux.cosphi;
  const double sinla = aux.sinla;
  const double cosla = aux.cosla;
  const double costwola = cosla * cosla - sinla * sinla;
  const double sintwola = 2e0 * cosla * sinla;

  // Compute the normalized position vector of the Moon.
  const double rmon2 = aux.rmon * aux.rmon;

  // Compute the normalized position vector of the Sun.
  const double rsun2 = aux.rsun * aux.rsun;

  // (minor modification) compute some helpfull intermediate quantities,
  // to reduce the following computation lines.
  const double xs0m1 = xsun.x() * xsun.x() - xsun.y() * xsun.y();
  const double xm0m1 = xmon.x() * xmon.x() - xmon.y() * xmon.y();
  const double fac2sun = aux.fac2sun;
  const double fac2mon = aux.fac2mon;

  const double drsun =
      -3e0 / 4e0 * dhi * cosphi * cosphi * fac2sun *
      (xs0m1 * sintwola - 2e0 * xsun.x() * xsun.y() * costwola) / rsun2;

  const double drmon =
      -3e0 / 4e0 * dhi * cosphi * cosphi * fac2mon *
      (xm0m1 * sintwola - 2e0 * xmon.x() * xmon.y() * costwola) / rmon2;

  const double dnsun =
      3e0 / 2e0 * dli * sinphi * cosphi * fac2sun *
      (xs0m1 * sintwola - 2e0 * xsun.x() * xsun.y() * costwola) / rsun2;

  const double dnmon =
      3e0 / 2e0 * dli * sinphi * cosphi * fac2mon *
      (xm0m1 * sintwola - 2e0 * xmon.x() * xmon.y() * costwola) / rmon2;

  const double desun =
      -3e0 / 2e0 * dli * cosphi * fac2sun *
      (xs0m1 * costwola + 2e0 * xsun.x() * xsun.y() * sintwola) / rsun2;

  const double demon =
      -3e0 / 2e0 * dli * cosphi * fac2mon *
      (xm0m1 * costwola + 2e0 * xmon.x() * xmon.y() * sintwola) / rmon2;

  const double dr = drsun + drmon;
  const double dn = dnsun + dnmon;
  const double de = desun + demon;

  // Compute the corrections for the station.
  return {{dr * cosla * cosphi - de * sinla - dn * sinphi * cosla,
           dr * sinla * cosphi + de * cosla - dn * sinphi * sinla,
           dr * sinphi + dn * cosphi}};

  // Finished
}

/// @details This function gives the out-of-phase corrections induced by
///          mantle anelasticity in the diurnal band.
///          This function is a translation/wrapper for the fortran ST1IDIU
///          subroutine, found here :
///          http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  xsta    Geocentric position of the station (Note 1)
/// @param[in]  xsun    Geocentric position of the Sun (Note 2)
/// @param[in]  xmon    Geocentric position of the Moon (Note 2)
/// @param[in]  fac2sun Degree 2 TGP factor for the Sun (Note 3)
/// @param[in]  fac2mon Degree 2 TGP factor for the Moon (Note 3)
/// @return xcorsta Out of phase station corrections for diurnal band
///                     (3d vector)
///
/// @note
///   -# The (IGS) station is in ITRF co-rotating frame. All coordinates are
///      expressed in meters, as arrays, i.e. [x,y,z].
///   -# The position is in Earth Centered Earth Fixed (ECEF) frame.  All
///      coordinates are expressed in meters, as arrays, i.e. [x,y,z].
///   -# The expressions are computed in the main program. TGP is the tide
///      generated potential. The units are inverse meters.
///   -# Status: Class 1
///   -# This fucnction is part of the package dehanttideinel, see
///      ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/
///
/// @version 31.07.2009
///
/// @cite iers2010,
///       Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, "Tidal station
///       displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
Vector3 st1idiu(/*const Vector3 &xsta,*/ const Vector3 &xsun,
                const Vector3 &xmon, const TideAux &aux) noexcept {
  constexpr const double dhi = -0.0025e0;
  constexpr const double dli = -0.0007e0;

  // Compute the normalized position vector of the station.
  // const double rsta = aux.rsta;
  const double sinphi = aux.sinphi;
  const double cosphi = aux.cosphi;
  const double cos2phi = cosphi * cosphi - sinphi * sinphi;
  const double sinla = aux.sinla;
  const double cosla = aux.cosla;

  // Compute the normalized position vector of the Moon.
  const double rmon2 = aux.rmon * aux.rmon;

  // Compute the normalized position vector of the Sun.
  const double rsun2 = aux.rsun * aux.rsun;
  const double fac2sun = aux.fac2sun;
  const double fac2mon = aux.fac2mon;

  const double drsun = -3e0 * dhi * sinphi * cosphi * fac2sun * xsun.z() *
                       (xsun.x() * sinla - xsun.y() * cosla) / rsun2;

  const double drmon = -3e0 * dhi * sinphi * cosphi * fac2mon * xmon.z() *
                       (xmon.x() * sinla - xmon.y() * cosla) / rmon2;

  const double dnsun = -3e0 * dli * cos2phi * fac2sun * xsun.z() *
                       (xsun.x() * sinla - xsun.y() * cosla) / rsun2;

  const double dnmon = -3e0 * dli * cos2phi * fac2mon * xmon.z() *
                       (xmon.x() * sinla - xmon.y() * cosla) / rmon2;

  const double desun = -3e0 * dli * sinphi * fac2sun * xsun.z() *
                       (xsun.x() * cosla + xsun.y() * sinla) / rsun2;

  const double demon = -3e0 * dli * sinphi * fac2mon * xmon.z() *
                       (xmon.x() * cosla + xmon.y() * sinla) / rmon2;

  const double dr = drsun + drmon;
  const double dn = dnsun + dnmon;
  const double de = desun + demon;

  // Compute the corrections for the station.
  return Vector3{{dr * cosla * cosphi - de * sinla - dn * sinphi * cosla,
                  dr * sinla * cosphi + de * cosla - dn * sinphi * sinla,
                  dr * sinphi + dn * cosphi}};

  // Finished
}

/// @details This function gives the in-phase and out-of-phase corrections
///          induced by mantle anelasticity in the diurnal band.
///          This function is a translation/wrapper for the fortran STEP2DIU
///          subroutine, found here at
///          http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  xsta    Geocentric position of the IGS station (Note 1)
/// @param[in]  fhr     Fractional hours in the day (Note 2)
/// @param[in]  t       Centuries since J2000
/// @return xcorsta In phase and out of phase station corrections
///                     for diurnal band (Note 4)
///
/// @note
///    -# The IGS station is in ITRF co-rotating frame. All coordinates are
///       expressed in meters, as arrays, i.e. [x,y,z].
///    -# The fractional hours in the day is computed as the hour + minutes/60.0
///       + sec/3600.0.  The unit is expressed in Universal Time (UT).
///    -# ----
///    -# All coordinates are expressed in meters.
///    -# Status: Class 1
///    -# This fucnction is part of the package dehanttideinel, see
///       ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/
///
/// @version 20.10.2010
///
/// @cite iers2010,
///       Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, "Tidal station
///       displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477,
Vector3 step2diu(const Vector3 &xsta, const Step2Angles &angles,
                 const TideAux &aux) noexcept {

  constexpr double datdi[][9] = {
      {-3e0, 0e0, 2e0, 0e0, 0e0, -0.01e0, 0e0, 0e0, 0e0},
      {-3e0, 2e0, 0e0, 0e0, 0e0, -0.01e0, 0e0, 0e0, 0e0},
      {-2e0, 0e0, 1e0, -1e0, 0e0, -0.02e0, 0e0, 0e0, 0e0},
      {-2e0, 0e0, 1e0, 0e0, 0e0, -0.08e0, 0e0, -0.01e0, 0.01e0},
      {-2e0, 2e0, -1e0, 0e0, 0e0, -0.02e0, 0e0, 0e0, 0e0},
      {-1e0, 0e0, 0e0, -1e0, 0e0, -0.10e0, 0e0, 0e0, 0e0},
      {-1e0, 0e0, 0e0, 0e0, 0e0, -0.51e0, 0e0, -0.02e0, 0.03e0},
      {-1e0, 2e0, 0e0, 0e0, 0e0, 0.01e0, 0e0, 0e0, 0e0},
      {0e0, -2e0, 1e0, 0e0, 0e0, 0.01e0, 0e0, 0e0, 0e0},
      {0e0, 0e0, -1e0, 0e0, 0e0, 0.02e0, 0e0, 0e0, 0e0},
      {0e0, 0e0, 1e0, 0e0, 0e0, 0.06e0, 0e0, 0e0, 0e0},
      {0e0, 0e0, 1e0, 1e0, 0e0, 0.01e0, 0e0, 0e0, 0e0},
      {0e0, 2e0, -1e0, 0e0, 0e0, 0.01e0, 0e0, 0e0, 0e0},
      {1e0, -3e0, 0e0, 0e0, 1e0, -0.06e0, 0e0, 0e0, 0e0},
      {1e0, -2e0, 0e0, -1e0, 0e0, 0.01e0, 0e0, 0e0, 0e0},
      {1e0, -2e0, 0e0, 0e0, 0e0, -1.23e0, -0.07e0, 0.06e0, 0.01e0},
      {1e0, -1e0, 0e0, 0e0, -1e0, 0.02e0, 0e0, 0e0, 0e0},
      {1e0, -1e0, 0e0, 0e0, 1e0, 0.04e0, 0e0, 0e0, 0e0},
      {1e0, 0e0, 0e0, -1e0, 0e0, -0.22e0, 0.01e0, 0.01e0, 0e0},
      {1e0, 0e0, 0e0, 0e0, 0e0, 12.00e0, -0.80e0, -0.67e0, -0.03e0},
      {1e0, 0e0, 0e0, 1e0, 0e0, 1.73e0, -0.12e0, -0.10e0, 0e0},
      {1e0, 0e0, 0e0, 2e0, 0e0, -0.04e0, 0e0, 0e0, 0e0},
      {1e0, 1e0, 0e0, 0e0, -1e0, -0.50e0, -0.01e0, 0.03e0, 0e0},
      {1e0, 1e0, 0e0, 0e0, 1e0, 0.01e0, 0e0, 0e0, 0e0},
      {0e0, 1e0, 0e0, 1e0, -1e0, -0.01e0, 0e0, 0e0, 0e0},
      {1e0, 2e0, -2e0, 0e0, 0e0, -0.01e0, 0e0, 0e0, 0e0},
      {1e0, 2e0, 0e0, 0e0, 0e0, -0.11e0, 0.01e0, 0.01e0, 0e0},
      {2e0, -2e0, 1e0, 0e0, 0e0, -0.01e0, 0e0, 0e0, 0e0},
      {2e0, 0e0, -1e0, 0e0, 0e0, -0.02e0, 0e0, 0e0, 0e0},
      {3e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
      {3e0, 0e0, 0e0, 1e0, 0e0, 0e0, 0e0, 0e0, 0e0}};

  // Compute the phase angles in degrees.
  const double s = angles.s;
  const double tau = angles.tau;
  const double h = angles.h;
  const double p = angles.p;
  const double zns = angles.zns;
  const double ps = angles.ps;

  // const double rsta = aux.rsta;
  const double sinphi = aux.sinphi;
  const double cosphi = aux.cosphi;
  const double cosla = aux.cosla;
  const double sinla = aux.sinla;
  const double zla = std::atan2(xsta.y(), xsta.x());

  Vector3 xcorsta{{0e0, 0e0, 0e0}};

  double thetaf, dr, dn, de;
  const double f1 = 2e0 * sinphi * cosphi;
  const double g1 = cosphi * cosphi - sinphi * sinphi;
  for (int j = 0; j < 31; j++) {
    // Convert from degrees to radians.
    thetaf = (tau + datdi[j][0] * s + datdi[j][1] * h + datdi[j][2] * p +
              datdi[j][3] * zns + datdi[j][4] * ps) *
             DEG2RAD;

    const double stz = std::sin(thetaf + zla);
    const double ctz = std::cos(thetaf + zla);
    dr = datdi[j][5] * f1 * stz + datdi[j][6] * f1 * ctz;

    dn = datdi[j][7] * g1 * stz + datdi[j][8] * g1 * ctz;

    // DE=DATDI(8,J)*SINPHI*COS(THETAF+ZLA)+
    // Modified 20 June 2007
    de = datdi[j][7] * sinphi * ctz - datdi[j][8] * sinphi * stz;

    xcorsta.data[0] += dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
    xcorsta.data[1] += dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
    xcorsta.data[2] += dr * sinphi + dn * cosphi;
  }

  // Finished
  return xcorsta * 1e-3;
}

/// @details This function gives the in-phase and out-of-phase corrections
///          induced by mantle anelasticity in the long period band.
///          This function is a translation/wrapper for the fortran STEP2LON
///          subroutine, found here :
///          http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  xsta    Geocentric position of the IGS station (Note 1)
/// @param[in]  t       Centuries since J2000
/// @param[out] xcorsta In phase and out of phase station corrections
///                     for diurnal band (Note 2)
///
/// @note
///    -# The IGS station is in ITRF co-rotating frame. All coordinates are
///       expressed in meters, as arrays, i.e. [x,y,z].
///    -# All coordinates are expressed in meters.
///    -# Status: Class 1
///    -# This fucnction is part of the package dehanttideinel, see
///       ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/
///
/// @version 20.10.2010
///
/// @cite iers2010,
///       Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, "Tidal station
///       displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477,
Vector3 step2lon(/*const Vector3 &xsta,*/ const Step2Angles &angles,
                 const TideAux &aux) noexcept {

  constexpr double datdi[][9] = {
      {0e0, 0e0, 0e0, 1e0, 0e0, 0.47e0, 0.23e0, 0.16e0, 0.07e0},
      {0e0, 2e0, 0e0, 0e0, 0e0, -0.20e0, -0.12e0, -0.11e0, -0.05e0},
      {1e0, 0e0, -1e0, 0e0, 0e0, -0.11e0, -0.08e0, -0.09e0, -0.04e0},
      {2e0, 0e0, 0e0, 0e0, 0e0, -0.13e0, -0.11e0, -0.15e0, -0.07e0},
      {2e0, 0e0, 0e0, 1e0, 0e0, -0.05e0, -0.05e0, -0.06e0, -0.03e0}};

  // Compute the phase angles in degrees.
  // Compute the phase angles in degrees.
  const double s = angles.s;
  // const double tau = angles.tau;
  const double h = angles.h;
  const double p = angles.p;
  const double zns = angles.zns;
  const double ps = angles.ps;

  // const double rsta = aux.rsta;
  const double sinphi = aux.sinphi;
  const double cosphi = aux.cosphi;
  const double cosla = aux.cosla;
  const double sinla = aux.sinla;

  // double dr_tot = 0e0, dn_tot = 0e0;
  Vector3 xcorsta{{0e0, 0e0, 0e0}};
  double thetaf, dr, dn, de;
  const double f1 = (3e0 * sinphi * sinphi - 1e0) / 2e0;
  const double f2 = cosphi * sinphi * 2e0;
  for (int j = 0; j < 5; j++) {

    thetaf = (datdi[j][0] * s + datdi[j][1] * h + datdi[j][2] * p +
              datdi[j][3] * zns + datdi[j][4] * ps) *
             DEG2RAD;

    const double ct = std::cos(thetaf);
    const double st = std::sin(thetaf);

    dr = datdi[j][5] * f1 * ct + datdi[j][7] * f1 * st;

    dn = datdi[j][6] * f2 * ct + datdi[j][8] * f2 * st;

    de = 0e0;

    // dr_tot += dr;
    // dn_tot += dn;

    xcorsta.data[0] += dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
    xcorsta.data[1] += dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
    xcorsta.data[2] += dr * sinphi + dn * cosphi;
  }

  // Finished.
  return xcorsta * 1e-3;
}

/// @details This function computes the station tidal displacement
///          caused by lunar and solar gravitational attraction (see
///          References). The computations are calculated by the following
///          steps:<br> <b>Step 1):</b> General degree 2 and degree 3
///          corrections
///          + CALL ST1IDIU + CALL ST1ISEM + CALL ST1L1.<br>
///          <b>Step 2):</b> CALL STEP2DIU + CALL STEP2LON<br>
///          It has been decided that the <b>Step 3</b> non-correction for
///          permanent tide would not be applied in order to avoid a jump in the
///          reference frame. This Step 3 must be added in order to get the
///          non-tidal station position and to conform with the IAG Resolution.
///          This function is a translation/wrapper for the fortran
///          DEHANTTIDEINEL subroutine, found here :
///          http://maia.usno.navy.mil/conv2010/software.html
///
/// @param[in]  xsta   Geocentric position of the station (Note 1)
/// @param[in]  xsun   Geocentric position of the Sun (Note 2)
/// @param[in]  xmon   Geocentric position of the Moon (Note 2)
/// @param[in]  yr     Year (Note 3)
/// @param[in]  month  Month (Note 3)
/// @param[in]  day    Day of Month (Note 3)
/// @param[in]  fhr    Hour in the day (Notes 3 and 4)
/// @return dxtide Displacement vector (Note 5)
/// @return            Always 0.
///
/// @note
///     -# The station is in ITRF co-rotating frame.  All coordinates,
///        X, Y, and Z, are expressed in meters.
///     -# The position is in Earth Centered Earth Fixed (ECEF) frame.  All
///        coordinates are expressed in meters.
///     -# The values are expressed in Coordinated Universal Time (UTC).
///     -# The fractional hours in the day is computed as the hour +
///        minutes/60.0 + sec/3600.0.
///     -# The displacement vector is in the geocentric ITRF.  All components
///        are expressed in meters.
///     -# Parameters jc1 and jc2 constitute the date as Julian Centuries in TT
///        time scale. The actual date is given by the addition jc1+jc2.
///        Either jc1 or jc2 can be set to zero.
///     -# Status: Class 1
///     -# This fucnction is part of the package dehanttideinel, see
///        ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/
///
/// @version 19.12.2016
///
/// @cite iers2010,
///
///     Groten, E., 2000, Geodesists Handbook 2000, Part 4,
///     http://www.gfy.ku.dk/~iag/HB2000/part4/groten.htm. See also
///     ''Parameters of Common Relevance of Astronomy, Geodesy, and
///     Geodynamics," J. Geod., 74, pp. 134-140
///
///     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
///     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
///
///     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
///     IERS Technical Note No. 36, BKG (2010)
///
///     Pitjeva, E. and Standish, E. M., 2009, ''Proposals for the masses
///     of the three largest asteroids, the Moon-Earth mass ratio and the
///     Astronomical Unit," Celest. Mech. Dyn. Astr., 103, pp. 365-372
///
///     Ries, J. C., Eanes, R. J., Shum, C. K. and Watkins, M. M., 1992,
///     ''Progress in the Determination of the Gravitational Coefficient
///     of the Earth," Geophys. Res. Lett., 19(6), pp. 529-531
Vector3
iers2010::dehanttideinel_impl(const Vector3 &xsta, const Vector3 &xsun,
                              const Vector3 &xmon,
                              dso::datetime<dso::milliseconds> epoch) noexcept {
  // nominal second degree and third degree love numbers and shida numbers
  constexpr const double h20 = 0.6078e0;
  constexpr const double l20 = 0.0847e0;
  constexpr const double h3 = 0.292e0;
  constexpr const double l3 = 0.015e0;

  // scalar product of station vector with sun/moon vector
  TideAux aux(&xsta, &xsun, &xmon, 0e0, 0e0);
  const double rmon = aux.rmon;
  const double rsta = aux.rsta;
  const double rsun = aux.rsun;
  const double scm = xsta.dot_product(xmon);
  const double scs = xsta.dot_product(xsun);
  const double scsun = scs / rsta / rsun;
  const double scmon = scm / rsta / rmon;

  // computation of new h2 and l2 (Equation 2, Chapter 7)
  const double cosphi = aux.cosphi;
  const double h2 = h20 - 0.0006e0 * (1e0 - 3e0 / 2e0 * cosphi * cosphi);
  const double l2 = l20 + 0.0002e0 * (1e0 - 3e0 / 2e0 * cosphi * cosphi);

  // P2 term (Equation 5, Chapter 7)
  const double p2sun = 3e0 * (h2 / 2e0 - l2) * scsun * scsun - h2 / 2e0;
  const double p2mon = 3e0 * (h2 / 2e0 - l2) * scmon * scmon - h2 / 2e0;

  // P3 term (Equation 6, Chapter 7)
  const double p3sun = 5e0 / 2e0 * (h3 - 3e0 * l3) * std::pow(scsun, 3) +
                       3e0 / 2e0 * (l3 - h3) * scsun;
  const double p3mon = 5e0 / 2e0 * (h3 - 3e0 * l3) * std::pow(scmon, 3) +
                       3e0 / 2e0 * (l3 - h3) * scmon;

  // term in direction of sun/moon vector
  const double x2sun = 3e0 * l2 * scsun;
  const double x2mon = 3e0 * l2 * scmon;
  const double x3sun = 3e0 * l3 / 2e0 * (5e0 * scsun * scsun - 1e0);
  const double x3mon = 3e0 * l3 / 2e0 * (5e0 * scmon * scmon - 1e0);

  // factors for sun/moon using iau current best estimates (see references)
  const double mass_ratio_sun = 332946.0482e0;
  const double mass_ratio_moon = 0.0123000371e0;
  const double re = /*6378136.6e0*/ iers2010::Re;
  const double fac2sun = mass_ratio_sun * re * std::pow(re / rsun, 3);
  const double fac2mon = mass_ratio_moon * re * std::pow(re / rmon, 3);
  const double fac3sun = fac2sun * (re / rsun);
  const double fac3mon = fac2mon * (re / rmon);

  // total displacement
  // for (int i = 0; i < 3; i++) {
  //  dxtide[i] = fac2sun * (x2sun * xsun[i] / rsun + p2sun * xsta[i] / rsta) +
  //              fac2mon * (x2mon * xmon[i] / rmon + p2mon * xsta[i] / rsta) +
  //              fac3sun * (x3sun * xsun[i] / rsun + p3sun * xsta[i] / rsta) +
  //              fac3mon * (x3mon * xmon[i] / rmon + p3mon * xsta[i] / rsta);
  //}
  Vector3 dxtide = fac2sun * (x2sun * xsun / rsun + p2sun * xsta / rsta) +
                   fac2mon * (x2mon * xmon / rmon + p2mon * xsta / rsta) +
                   fac3sun * (x3sun * xsun / rsun + p3sun * xsta / rsta) +
                   fac3mon * (x3mon * xmon / rmon + p3mon * xsta / rsta);
  aux.fac2mon = fac2mon;
  aux.fac2sun = fac2sun;

  // corrections for the out-of-phase part of love numbers (part h_2^(0)i
  // and l_2^(0)i )
  // --------------------------------------------------------------------

  // first, for the diurnal band
  dxtide += st1idiu(/*xsta, */ xsun, xmon, aux);

  // second, for the semi-diurnal band
  dxtide += st1isem(/*xsta, */ xsun, xmon, aux);

  // corrections for the latitude dependence of love numbers (part l^(1) )
  dxtide += st1l1(/*xsta, */ xsun, xmon, aux);

  // CONSIDER CORRECTIONS FOR STEP 2
  // -------------------------------------------------------------------------

  // corrections for the diurnal band:
  // first, we need to know the date converted in julian centuries
  // UTC to TT time
  double fhr = epoch.sec().to_fractional_seconds();
  fhr /= 3600e0;
  int dat = dso::dat(epoch.mjd());
  epoch.add_seconds(dso::milliseconds(dat * 1e3) + dso::milliseconds(32184));
  double t = (epoch.as_mjd() + dso::mjd0_jd - dso::j2000_jd) / 36525e0;
  Step2Angles angles(t, fhr);

  //  second, we can call the subroutine step2diu, for the diurnal band
  //+ corrections, (in-phase and out-of-phase frequency dependence):
  dxtide += step2diu(xsta, angles, aux);

  //  corrections for the long-period band,
  //+ (in-phase and out-of-phase frequency dependence):
  dxtide += step2lon(/*xsta, */ angles, aux);

  // CONSIDER CORRECTIONS FOR STEP 3
  // --------------------------------------------------------------------------

  // uncorrect for the permanent tide
  /*
   * double sinphi = xsta[2]/rsta;
   * double cosphi = sqrt( xsta[0]*xsta[0]+xsta[1]*xsta[1] ) / rsta;
   * double cosla  = xsta[0]/cosphi/rsta;
   * double sinla  = xsta[1]/cosphi/rsta;
   * double dr     =
   * -sqrt(5e0/4e0/PI)*h2*0.31460e0*(3e0/2e0*sinphi*sinphi-0.5e0); double dn =
   * -sqrt(5e0/4e0/PI)*l2*0.31460e0*3e0*cosphi*sinphi; dxtide[0]    +=
   * -dr*cosla*cosphi+dn*cosla*sinphi; dxtide[1]    +=
   * -dr*sinla*cosphi+dn*sinla*sinphi; dxtide[2]    += -dr*sinphi -dn*cosphi;
   */

  // Finished.
  return dxtide;
}