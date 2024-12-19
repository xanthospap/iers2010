# Fundamental Arguments And Planetary Angles

In the concept of IERS standards, the fundamental arguments, or namely the Delaunay 
variables (see [IERS](#IERS2010)(Section 5.7)), are:

  * $l$  : Mean anomaly of the Moon [rad]
  * $l'$ : Mean anomaly of the Sun [rad]
  * $F$  : L - OM [rad]
  * $D$  : Mean elongation of the Moon from the Sun [rad]
  * $\Omega$ : Mean longitude of the ascending node of the Moon [rad]

Within the [IERS](#IERS2010), these are oftentime labelled $F_j$. The 
implementation for the computation of the fundamental arguments adopted here, 
closesly follows the [SOFA](#SOFA) library and is consistent with its results 
to at least $1e-9 [arcsec]$.

The $F_j$ are functions of time, and the **angular frequency** of the respective 
term is given by $\omega = d(F_j) / dt$

## Testing

## Issues

When the compiler optimizations are turned on (i.e. for the production build), the 
implementation of the fundamental arguments produces way larger discrepancies w.r.t 
the SOFA library (i.e. in the order of $1e-1 [arcsec]$. For this reason, this 
part of the software is included withi (pre-proccessor) `PRAGMA` guards (see the 
file [fundarg.hpp](../blobl/cleanup/src/fundarg.hpp) and 
[fundarg.cpp](../blobl/cleanup/src/fundarg/fundarg.cpp). We should check whay 
this happens.

## References

<a name="IERS2010"></a>IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). 
(IERS Technical Note ; 36) Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie, 
2010. 179 pp., ISBN 3-89888-989-6

<a name="SOFA"></a>Software Routines from the IAU SOFA Collection were used. 
Copyright © International Astronomical Union Standards of Fundamental Astronomy
(http://www.iausofa.org)
