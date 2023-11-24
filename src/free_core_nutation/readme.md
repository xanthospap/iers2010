# Free Core Nutation

Free core nutation is a free retrograde diurnal motion of the Earth’s rotation axis with respect
to the Earth, which appears as a motion of the CIP in the GCRS. FCN is not included in 
the IAU 2000A nutation model and hence has to be taken into account depending on 
accuracy requirements. 

The model syggested by the [IERS](#IERS2010) is the The FCN model of [Lambert 2007](#Lambert-2007); 
this uses regularly updated coefficients that can be found at 
[Paris Observatory Geodetic VLBI Center](http://ivsopar.obspm.fr/fcn/).

Note that the unmodeled FCN motion of the CIP is included in the published IERS celestial pole
offsets dX and dY. **These offsets should not be applied when the FCN model is used.**

## Input File(s)

To use the Lambert model, we must have the (regularly) updated model coefficients. 
These are published (in ascii) at [Paris Observatory Geodetic VLBI Center](http://ivsopar.obspm.fr/fcn/).
This is handled by the

# References

<a name="Lambert-2007"></a>Lambert, S., 2007, “Empirical modeling of the retrograde Free Core Nutation,” available at
ftp://hpiers.obspm.fr/eop-pc/models/fcn/notice.pdf.

<a name="IERS2010"</a>IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). 
(IERS Technical Note ; 36) Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie, 
2010. 179 pp., ISBN 3-89888-989-6, Section 5.5.5

for the most stringent accuracy applications, a FCN model may be incorporated to
account for the FCN contribution to the CIP motion in the GCRS.
