#include "hardisp.hpp"
#include "iers2010.hpp"

/// @details This subroutine returns the frequency and phase of a tidal
///          constituent when its Doodson number is given as input.
///
/// @param[in]  idood Doodson number of a tidal constituent (6-element integer
/// array)
/// @param[in]  itm   Date as integer array in UTC. The format should be:
///                   [year, day_of_year, hours, minutes, seconds].
/// @param[out] freq  Frequency of a tidal constituent
/// @param[out] phase Phase of a tidal constituent (Note 1)
/// @return           Always 0.
///
/// @note
///     -# The phases must be decreased by 90 degrees if the sum of the order
///        and the species number is odd (as for the 2nd degree diurnals, and
///        3rd degree low frequency and semidiurnals).
///        These phases may need further adjustment to allow for the spherical
///        harmonic normalization used; e.g. for that used for the potential
///        by Cartwright and Tayler, 180 degrees must be added for (species,
///        order) = (1,2), (1,3), or (3,3).
///     -# Status:  Class 1 model
///
/// @version 11.09.2013
///
/// @cite iers2010
///
/// @todo There is a bug here! The line:
/// 'double delta  ( iers2010::hisp::etutc (year) );'
/// returns a delta != 0. The FORTRAN routine however gives delta=0.
/// WTF ??
int iers2010::hisp::tdfrph(const int idood[6],
                           ngpt::datetime<ngpt::seconds> iepoch, double &freq,
                           double &phase) {
  static ngpt::datetime<ngpt::seconds> last_epoch;
  static double d[6];
  static double dd[6];

  //  Test to see if time has changed; if so, set the phases and frequencies
  //  for each of the Doodson arguments
  if (iepoch != last_epoch) {

    ngpt::datetime<ngpt::milliseconds> epoch =
        iepoch.cast_to<ngpt::milliseconds>();
    // Convert times to Julian days (UT) then to Julian centuries
    // from J2000.0 (TT)
    int dat = ngpt::dat(epoch.mjd());
    double dayfr = epoch.sec().fractional_days();
    epoch.add_seconds(ngpt::milliseconds(dat * 1e3) +
                      ngpt::milliseconds(32184));
    double t = (epoch.as_mjd() + ngpt::mjd0_jd - ngpt::j2000_jd) / 36525e0;
    // printf("\ndat=%5d, dayfr=%15.10f, jd=%20.10f t=%15.10f", dat, dayfr,
    // epoch.as_mjd() + ngpt::mjd0_jd, t); t = 0.094799449636217e0; dayfr =
    // 0.049131944444444e0;

    // IERS expressions for the Delaunay arguments, in degrees
    double f1{134.9634025100e0 +
              t * (477198.8675605000e0 +
                   t * (0.0088553333e0 +
                        t * (0.0000143431e0 + t * (-0.0000000680e0))))};

    double f2{357.5291091806e0 +
              t * (35999.0502911389e0 +
                   t * (-0.0001536667e0 +
                        t * (0.0000000378e0 + t * (-0.0000000032e0))))};

    double f3{93.2720906200e0 +
              t * (483202.0174577222e0 +
                   t * (-0.0035420000e0 +
                        t * (-0.0000002881e0 + t * (0.0000000012e0))))};

    double f4{297.8501954694e0 +
              t * (445267.1114469445e0 +
                   t * (-0.0017696111e0 +
                        t * (0.0000018314e0 + t * (-0.0000000088e0))))};

    double f5{125.0445550100e0 +
              t * (-1934.1362619722e0 +
                   t * (0.0020756111e0 +
                        t * (0.0000021394e0 + t * (-0.0000000165e0))))};

    // Convert to Doodson (Darwin) variables
    d[0] = 360e0 * dayfr - f4;
    d[1] = f3 + f5;
    d[2] = d[1] - f4;
    d[3] = d[1] - f1;
    d[4] = -f5;
    d[5] = d[2] - f2;

    //  Find frequencies of Delauney variables (in cycles/day), and from
    //+ these the same for the Doodson arguments
    double fd1{0.0362916471e0 + 0.0000000013e0 * t};
    double fd2{0.0027377786e0};
    double fd3{0.0367481951e0 - 0.0000000005e0 * t};
    double fd4{0.0338631920e0 - 0.0000000003e0 * t};
    double fd5{-0.0001470938e0 + 0.0000000003e0 * t};
    dd[0] = 1e0 - fd4;
    dd[1] = fd3 + fd5;
    dd[2] = dd[1] - fd4;
    dd[3] = dd[1] - fd1;
    dd[4] = -fd5;
    dd[5] = dd[2] - fd2;

    // copy the just used date to the static one for next use of function
    last_epoch = iepoch;

  } // End of intialization (likely to be called only once)

  //  Compute phase and frequency of the given tidal constituent
  freq = 0e0;
  phase = 0e0;
  for (int i = 0; i < 6; i++) {
    freq += ((double)idood[i]) * dd[i];
    phase += ((double)idood[i]) * d[i];
  }

  // Adjust phases so that they fall in the positive range 0 to 360
  phase = std::fmod(phase, 360e0);
  if (phase < 0e0)
    phase += 360e0;

  // Finished
  return 0;
}
