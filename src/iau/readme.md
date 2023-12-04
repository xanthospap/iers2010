# Earth Rotation Angle (ERA)

Computation of ERA follows Equation 5.15 in [IERS Conventions 2010](#IERS2010). 
Differences between this implementation and SOFA, are smaller than $1e-9 [arcsec]$.
For testing, see [era00.cpp](../blob/cleanup/test/sofa/era00.cpp).


## References

<a name="IERS2010"></a>IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). 
(IERS Technical Note ; 36) Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie, 
2010. 179 pp., ISBN 3-89888-989-6
