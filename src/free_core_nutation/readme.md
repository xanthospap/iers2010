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
Parsing is handled by the `parse_lambert_coefficients` function. Once the model 
coefficients are parsed, a call to `lambert_fcn` can be used (with the 
coefficients and time as input parameters), to compute CIP corrections.

## IERS 2010 FCNNUT distributed routine

The IERS provides the FCNNUT.f file/subroutine to model FCN, using the Lambert 
model. It provides the needed coefficients within the source code. The most 
recent coefficients are up to 2013 (at the time of writing). To use this 
version of the model, you do not need any input files, just call 
`lambert_fcn` with time as input parameter.

## Tests

See the file [lambert.cpp](test/unit_tests/lambert.cpp); it checks results of 
the implemented model against the ones obtained from the routine the 
[IERS distributed FORTRAN routine](http://ivsopar.obspm.fr/fcn/fcndrv.f).

<div style="text-align: right"> Last update Dec-2023 </div>

## Issues

Well, quite frankly i am not sure what time scale should the time arguments be in. 
I **think** it is UTC, but i am not sure. For this reason, and probably because it 
pretty much does not make any difference, date/times are treated as continuous, naive 
epochs (hence non-UTC). At some point though, this should be clarified.

When clarified, both the instances of `fcn::LambertCoefficients` should be treated 
accordingly and the function parameters.

## References

<a name="Lambert-2007"></a>Lambert, S., 2007, “Empirical modeling of the retrograde Free Core Nutation,” available at
ftp://hpiers.obspm.fr/eop-pc/models/fcn/notice.pdf.

<a name="IERS2010"></a>IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). 
(IERS Technical Note ; 36) Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie, 
2010. 179 pp., ISBN 3-89888-989-6, Section 5.5.5
