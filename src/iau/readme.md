# Earth Rotation Angle (ERA)

Computation of ERA follows Equation 5.15 in [IERS Conventions 2010](#IERS2010). 
Differences between this implementation and SOFA, are smaller than $1e-9 [arcsec]$.
For testing, see [era00.cpp](../blob/cleanup/test/sofa/era00.cpp).

# TIO Locator $s\prime$

The quantity $s\prime$ (i.e. the TIO[^1] locator) represents the motion of the TIO 
in the ITRS. It is computed using the current mean amplitudes for the Chandlerian 
and annual wobbles [Lambert and Bizouard, 2002](#TIORefPaper), 
as suggested by [IERS Conventions 2010](#IERS2010) (Section 5.5.2):
$s\prime = -47 µas t$. 

[^1]: TIO was labelled TEO before IAU 2006 Resolution B2, which harmonized 
“intermediate” in the names of the pole and the origin (i.e. celestial and
terrestrial intermediate origins, CIO and TIO instead of CEO and TEO, respectively).

# Coordinates of CIP in the GCRS

The coordinates of the CIP in the GCRS (i.e. $(X,Y)$), are computed as developments 
as function of time (as described in [IERS Conventions 2010](#IERS2010) (Section 5.5.4)).
The developments are valid at the microarcsecond level, based on the IAU 2006 
precession and IAU 2000A nutation.

Differences between this implementation and SOFA (see function `iauXy06`), are smaller 
than $1e-9 [arcsec]$. For testing, see [xy06a.cpp](../blob/cleanup/test/sofa/xy06a.cpp).

VLBI observations have shown that there are deﬁciencies in the IAU 2006/2000A 
precession-nutation model of the order of $0.2 mas$, mainly due to the fact 
that the free core nutation (FCN) is not part of the model. The IERS publishes 
observed estimates of the corrections to the IAU precession-nutation model (w.r.t 
the conventional celestial pole position deﬁned by the models), named "celestial pole oﬀsets". 
Such time-dependent oﬀsets from the direction of the pole of the GCRS must be 
provided as corrections $\delta X$ and $\delta Y$ to the $X$ and $Y$ coordinates.
Using these oﬀsets, the corrected celestial position of the CIP is given by:
$X = X(IAU 2006/2000) + \delta X$ and $Y = Y(IAU 2006/2000) + \delta Y$.

## References

<a name="IERS2010"></a>IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). 
(IERS Technical Note ; 36) Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie, 
2010. 179 pp., ISBN 3-89888-989-6

<a name="TIORefPaper"></a>Lambert, S. and Bizouard, C., 2002, 
“Positioning the Terrestrial Ephemeris Origin in the International Terrestrial Reference Frame,” 
Astron. Astrophys., 394(1), pp. 317–321, doi:10.1051/0004-6361:20021139.
