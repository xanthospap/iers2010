#include "iau.hpp"

void iers2010::sofa::p06e(const dso::TwoPartDate &mjd_tt, double &eps0,
                          double &psia, double &oma, double &bpa, double &bqa,
                          double &pia, double &bpia, double &epsa, double &chia,
                          double &za, double &zetaa, double &thetaa, double &pa,
                          double &gam, double &phi, double &psi) noexcept {

  /* Interval between fundamental date J2000.0 and given date (JC). */
  const double t = mjd_tt.jcenturies_sinceJ2000();

  /* Obliquity at J2000.0. */
  eps0 = dso::sec2rad(84381.406e0);

  /* Luni-solar precession. */
  psia = dso::sec2rad(
      (5038.481507e0 +
       (-1.0790069e0 +
        (-0.00114045e0 + (0.000132851e0 + (-0.0000000951e0) * t) * t) * t) *
           t) *
      t);

  /* Inclination of mean equator with respect to the J2000.0 ecliptic. */
  oma =
      eps0 +
      dso::sec2rad(
          (-0.025754e0 +
           (0.0512623e0 +
            (-0.00772503e0 + (-0.000000467e0 + (0.0000003337e0) * t) * t) * t) *
               t) *
          t);

  /* Ecliptic pole x, J2000.0 ecliptic triad. */
  bpa = dso::sec2rad(
      (4.199094e0 +
       (0.1939873e0 +
        (-0.00022466e0 + (-0.000000912e0 + (0.0000000120e0) * t) * t) * t) *
           t) *
      t);

  /* Ecliptic pole -y, J2000.0 ecliptic triad. */
  bqa = dso::sec2rad(
      (-46.811015e0 +
       (0.0510283e0 +
        (0.00052413e0 + (-0.000000646e0 + (-0.0000000172e0) * t) * t) * t) *
           t) *
      t);

  /* Angle between moving and J2000.0 ecliptics. */
  pia = dso::sec2rad(
      (46.998973e0 +
       (-0.0334926e0 +
        (-0.00012559e0 + (0.000000113e0 + (-0.0000000022e0) * t) * t) * t) *
           t) *
      t);

  /* Longitude of ascending node of the moving ecliptic. */
  bpia = dso::sec2rad(
      629546.7936e0 +
      (-867.95758e0 +
       (0.157992e0 +
        (-0.0005371e0 + (-0.00004797e0 + (0.000000072e0) * t) * t) * t) *
           t) *
          t);

  /* Mean obliquity of the ecliptic. */
  epsa = iers2010::sofa::obl06(mjd_tt);

  /* Planetary precession. */
  chia = dso::sec2rad(
      (10.556403e0 +
       (-2.3814292e0 +
        (-0.00121197e0 + (0.000170663e0 + (-0.0000000560e0) * t) * t) * t) *
           t) *
      t);

  /* Equatorial precession: minus the third of the 323 Euler angles. */
  za = dso::sec2rad(
      -2.650545e0 +
      (2306.077181e0 +
       (1.0927348e0 +
        (0.01826837e0 + (-0.000028596e0 + (-0.0000002904e0) * t) * t) * t) *
           t) *
          t);

  /* Equatorial precession: minus the first of the 323 Euler angles. */
  zetaa = dso::sec2rad(
      (2.650545e0 +
       (2306.083227e0 +
        (0.2988499e0 +
         (0.01801828e0 + (-0.000005971e0 + (-0.0000003173e0) * t) * t) * t) *
            t) *
           t));

  /* Equatorial precession: second of the 323 Euler angles. */
  thetaa = dso::sec2rad(
      (2004.191903e0 +
       (-0.4294934e0 +
        (-0.04182264e0 + (-0.000007089e0 + (-0.0000001274e0) * t) * t) * t) *
           t) *
      t);

  /* General precession. */
  pa = dso::sec2rad(
      (5028.796195e0 +
       (1.1054348e0 +
        (0.00007964e0 + (-0.000023857e0 + (-0.0000000383e0) * t) * t) * t) *
           t) *
      t);

  /* Fukushima-Williams angles for precession. */
  gam = dso::sec2rad(
      (10.556403e0 +
       (0.4932044e0 +
        (-0.00031238e0 + (-0.000002788e0 + (0.0000000260e0) * t) * t) * t) *
           t) *
      t);

  phi =
      eps0 +
      dso::sec2rad(
          (-46.811015e0 +
           (0.0511269e0 +
            (0.00053289e0 + (-0.000000440e0 + (-0.0000000176e0) * t) * t) * t) *
               t) *
          t);

  psi = dso::sec2rad(
      (5038.481507e0 +
       (1.5584176e0 +
        (-0.00018522e0 + (-0.000026452e0 + (-0.0000000148e0) * t) * t) * t) *
           t) *
      t);

  /* Finished. */
  return;
}
