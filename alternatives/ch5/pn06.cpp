#include "iau.hpp"
#include "iersc.hpp"

void iers2010::sofa::pn06(double date1, double date2, double dpsi, double deps,
                          double &epsa, iers2010::RotationMatrix3 &rb,
                          iers2010::RotationMatrix3 &rp,
                          iers2010::RotationMatrix3 &rbp,
                          iers2010::RotationMatrix3 &rn,
                          iers2010::RotationMatrix3 &rbpn) noexcept {
  // Bias-precession Fukushima-Williams angles of J2000.0 = frame bias. 
  double gamb, phib, psib, eps;
  iers2010::sofa::pfw06(iers2010::DJM0, iers2010::DJM00, gamb, phib, psib, eps);

  // B matrix. 
  rb = iers2010::sofa::fw2m(gamb, phib, psib, eps);

  // Bias-precession Fukushima-Williams angles of date. 
  iers2010::sofa::pfw06(date1, date2, gamb, phib, psib, eps);

  // Bias-precession matrix. 
  rbp = iers2010::sofa::fw2m(gamb, phib, psib, eps);

  // Solve for precession matrix.
  rp = rbp * rb.transpose();

  // Equinox-based bias-precession-nutation matrix. 
  rbpn = iers2010::sofa::fw2m(gamb, phib, psib + dpsi, eps + deps);

  // Solve for nutation matrix. 
  rn = rbpn * rbp.transpose();

  // Obliquity, mean of date. 
  epsa = eps;

  // Finished.
  return;
}