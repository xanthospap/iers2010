# Transformation between the International Terrestrial Reference System and the Geocentric Celestial Reference System

The transformation to be used to relate the ITRS to the GCRS at the date $t$ 
of the observation can be written as:

$$ \mathbf{x}_{GCRS} = Q(t) R(t) W(t) \mathbf{x}_{ITRS} $$

where:
 - $Q(t)$ accounts for *precession* and *nutation*,
 - $R(t)$ accounts for Earth's rotation, and
 - $W(t)$ accounts for polar motion

## Precession, Nutation and Matrix $Q(t)$

The matrix $Q(t)$ accounts for *precession* and *nutation* and represents the 
dynamical motion of the **Celestial Intermediate Pole** in response to torques 
arising from the gravitational attraction of the Sun, Moon and planets (see 
[Almanac](#EXAlmanac), Section 6.1). Traditionally, precession accounted for the 
secular part of the motion, while nutation for the periodic component and their 
impact was modeled using two individual (rotation) matrices $P(t)$ and 
$N(t)$.

The matrix $Q(t)$, relates (i.e. transforms between) the CIRS and GCRS frames.
It can be computed by the formulae (see [IERS Conventions 2010](#IERS2010) 
(Eq. 5.6 through 5.10)):

$$ $Q(t) = R_3 (-E) R_2 (-d) R_3 (E) R_3 (s) $$,

or

$$ Q(t) = \begin{pmatrix} 1-\alpha X^2 & -\alpha XY   & X \\
                                -\alpha XY  & 1-\alpha Y^2 & Y \\
                                -X          & -Y           & 1-\alpha (X^2 + Y^2)
               \end{pmatrix} \cdot R_3 (s) $$

where:
    - $X = \sin d \cos E$,
    - $Y = \sin d \sin E$,
    - $Z = \cos d$, 
    - $s$ is the ``CIO locator'', and
    - $\alpha = 1/(1 + \cos d ) \approx 1/2 + 1/8(X^2 + Y^2)$ within $1 [\mu as]$

$X$ and $Y$ are the [CIP coordinates](#cipingcrs).

## Coordinates of CIP in the GCRS <a name="cipingcrs"></a>

The coordinates of the CIP in the GCRS (i.e. $(X,Y)$), are computed as developments 
as function of time (as described in [IERS Conventions 2010](#IERS2010) (Section 5.5.4)).
The developments are valid at the microarcsecond level, based on the IAU 2006 
precession and IAU 2000A nutation.

Differences between this implementation and SOFA (see function `iauXy06`), are 
way below $1e-12 [arcsec]$ (i.e. machine prescision). For testing, see 
[xy06a.cpp](../blob/cleanup/test/sofa/xy06a.cpp).

VLBI observations have shown that there are deﬁciencies in the IAU 2006/2000A 
precession-nutation model of the order of $0.2 mas$, mainly due to the fact 
that the free core nutation (FCN) is not part of the model. The IERS publishes 
observed estimates of the corrections to the IAU precession-nutation model (w.r.t 
the conventional celestial pole position deﬁned by the models), named "celestial pole oﬀsets". 
Such time-dependent oﬀsets from the direction of the pole of the GCRS must be 
provided as corrections $\delta X$ and $\delta Y$ to the $X$ and $Y$ coordinates.
Using these oﬀsets, the corrected celestial position of the CIP is given by:
$X = X(IAU 2006/2000) + \delta X$ and $Y = Y(IAU 2006/2000) + \delta Y$.


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

# CIO Locator $s$

For the computation of the CIO position, $s$, we use the IAU 2006/2000A precession-nutation 
model (see [IERS Conventions 2010](#IERS2010) (Section 5.5.6)). Note that therein, 
two approaches are outlined, i.e based on Table 5.2d and one on Table 5.2c, the 
latter being of slightly better accuracy. In this implementation we use an implementation 
based on Table 5.2c, i.e. a development of $s(t) + X_{CIP}Y_{CIP} /2$ including all 
terms larger than $0.1 \mu as$.

Differences between this implementation and SOFA, are smaller than $1e-9 [arcsec]$.
For testing, see [s06.cpp](../blob/cleanup/test/sofa/s06.cpp).

Due to the fact that the need to compute the CIO locator $s$ is coupled with the 
CIP coordinates $(X,Y)$, and both computations require linusolar and planetary 
arguments, the latter can be directly included in the computation by a dedicated 
parameter. E.g.
```
using namespace dso;

MjdEpoch tt(....);
#
/* Store lunisolar and planetary args for epoch tt in lpargs array. */
double lpargs[14];

/* Compute X and Y CIP coordinates in [rad], and store lunisolar and planetary 
 * arguments in lpargs.
 */
xycip06a(tt, xcip, ycip, lpargs);

/* Compute the CIO locator, s; do not recompute lunisolar and planetary args. 
 * Note that the computations are performed for the **same** epoch, tt.
 */
double s = s06(tt, xcip, ycip, lpargs);
```

## References

<a name="IERS2010"></a>IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). 
(IERS Technical Note ; 36) Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie, 
2010. 179 pp., ISBN 3-89888-989-6

<a name="TIORefPaper"></a>Lambert, S. and Bizouard, C., 2002, 
“Positioning the Terrestrial Ephemeris Origin in the International Terrestrial Reference Frame,” 
Astron. Astrophys., 394(1), pp. 317–321, doi:10.1051/0004-6361:20021139.

<a name="EXAlmanac"></a>A Seidelmann, P. Kenneth, Urban, Sean E., 
Explanatory Supplement to the Astronomical Almanac, Third Edition,
American Astronomical Society Meeting Abstracts, 2010
