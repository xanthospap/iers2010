## GCRS to ITRS

The information given below, are quoted from [^1], Section 5.9 and are relevant 
to the adoption/translation of the corresponding SOFA routines within the libiers2010 
library.

[^1]: Gérard Petit and Brian J. Luzum, editors. IERS Conventions (2010), volume 36 of
IERS Technical Note, 2010. International Earth Rotation and Reference Systems
Service (IERS), International Earth Rotation and Reference Systems Service
(IERS).

| Routine       | Brief Description                                           | IAU Resolution 
| :------------ | :----------------                                           | :-------------
| c2i06a        | CIO: celestial to intermediate matrix                       | IAU 2006/2000A (Q^-1)
| c2t06a        | CIO: celestial to terrestrial matrix                        | IAU 2006/2000A
| pnm06a        | Equinox: celestial to true matrix                           | IAU 2006/2000A (PNB)

| Routine       | Brief Description                                           | IAU Resolution | Equation from [^1] | Test Against SOFA/IAU
| :------------ | :----------------                                           | :------------- | :----------------- | :--------------------
| ee06a         | Equation of the equinoxes                                   | IAU 2006/2000A |                    | [x]
| era00         | Earth Rotation Angle                                        | IAU 2000       | 5.15               | [x]
| gst06a        | Greenwich (apparent) Sidereal Time                          | IAU 2006/2000A |                    | [x]
| gmst06a       | Greenwich Mean Sidereal Time                                | IAU 2006/2000A |                    | [x]
| nut06a        | Nutation components                                         | IAU 2006/2000A |                    | [x]
| pfw06         | 4-angle Fukushima-Williams precession angles                | IAU 2006       |                    | [x]
| pn06          | Bias, precession, nutation matrices given Δψ, Δε (B, P, N, NPB) | IAU 2006/2000A |                | [ ]
| pom00         | Polar motion matrix                                         | IAU 2006/2000A |                    | [x]
| s06           | CIO locator s, given X and Y                                | IAU 2006       |                    | [x]
| s06a          | CIO locator s                                               | IAU 2006/2000A |                    | [-]
| sp00          | TIO locator s'                                              | IAU 2000       |                    | [ ]
| xy06          | X, Y from semi-analytical series                            | IAU 2006/2000A |                    | [ ]
| xys06a        | X, Y, s                                                     | IAU 2006/2000A |                    | [ ]

| Routine       | Brief Description                                                 |
| :------------ | :----------------                                                 |
| c2ixys        | CIO: celestial to intermediate matrix, given X, Y and s (Q^-1)    |
| c2tcio        | CIO: celestial to terrestrial matrix                              |
| c2teqx        | Equinox: celestial to terrestrial matrix                          |
| fw2m          | Fukushima-Williams angles to rotation matrix (B, PB, or NPB)      |

The list above is by no means an exhaustive one. Nor is it the only solution to transform between
ITRS and GCRS. The matrix for the combined effects of nutation, precession and frame bias is Q(t)
in expression (5.1). For the CIO based transformation, this is the intermediate-to-celestial matrix,
it can be obtained (as the transpose) using the routine `C2IXYS`, starting from the CIP
position X,Y and the quantity s that defines the position of the CIO. The IAU 2006/2000A X, Y, s 
are available by calling the routine `XYS06A`. In the case of the equinox based transformation,
the counterpart to matrix Q(t) is the true-to-celestial matrix. To obtain this matrix requires
the nutation components ∆ψ and ∆ε; these can be predicted using the IAU 2000A model, with
adjustments to match IAU 2006 precession, by means of the routine `NUT06A`. Faster, but
less accurate, predictions are available from the `NUT00B` routine, which implements the IAU 2000B
truncated model. Once ∆ψ and ∆ε are known, the true-to-celestial matrix can be obtained by
calling the routine `PN06` and taking the transpose.

The intermediate component is the angle for Earth rotation that defines matrix R(t) in expression
(5.1). For the CIO based transformation, the angle in question is the Earth Rotation Angle, ERA,
which can be obtained by calling the routine `ERA00`. The counterpart in the case of the
equinox based transformation is the Greenwich (apparent) Sidereal Time. This can be obtained
by calling the routine `GST06`, given the celestial-to-true matrix that was obtained earlier.

The three components - the precession-nutation matrix, the Earth rotation quantity and the polar
motion matrix – are then assembled into the final terrestrial-to-celestial matrix by means of the
routine `C2TCIO` (CIO based) or `C2TEQX` (equinox based), and getting thier transposes as required.

Two methods to generate the terrestrial-to-celestial (i.e. ITRS-to-GCRS) matrix Q(t)R(t)W(t),
given TT and UT1, are set out below. In each case it is assumed that observed small corrections
to the IAU 2006/2000A model, either as ∆X, ∆Y or as d∆ψ, d∆ε, are available and need to be
included.

## Method (1): the CIO based transformation

The CIO based transformation is a function of the CIP coordinates X, Y and the quantity s.
For the given TT, call the routine `XY06` to obtain the IAU 2006/2000A X, Y from series
(see Section 5.5.4) and then the routine `S06` to obtain s. Any CIP corrections ∆X, ∆Y can now
be applied, and the corrected X, Y, s can be used to call the routine `C2IXYS`, giving the GCRS-to-
CIRS matrix. Next call the routine `ERA00` to obtain the ERA corresponding to the current UT1,
and apply it as an R3 rotation using the routine RZ, to form the CIRS-to-TIRS matrix. Given
xp , yp , and obtaining s0 by calling the routine `SP00`, the polar motion matrix (i.e. TIRS-to-ITRS)
is then produced by the routine `POM00`. The product of the two matrices (GCRS-to-TIRS and
TIRS-to-ITRS), (obtained by calling the routine `RXR` [^2]), is the GCRS-to-ITRS matrix (which can be
inverted by calling the routine `TR` to give the final result [^2]).

[^2]: The routines `RXR` and `TR` are SOFA routines not implemented here. Instead, the 
libiers2010 approach is to define a (kinda) matrix class (`dso::Mat3x3`) and return 
such matrices (instead of multi-dimensional arrays). Matrix multiplication and transpotition 
are defined for the instances of tis class.

## Method (2): the equinox based transformation

The classical transformation, based on angles and using sidereal time is also available.
Given TT, the IAU 2006/2000A nutation components ∆ψ, ∆ε are obtained by calling the 
routine `NUT06A`. Any corrections d∆ψ, d∆ε can now be applied. Next, the GCRS-to-true matrix
is obtained using the routine `PN06` (which employs the 4-rotation Fukushima-Williams method
described in Section 5.3.4, final paragraph). The classical GCRS-to-true matrix can also be gener-
ated by combining separate frame bias, precession and nutation matrices. The routine `BI00`
can be called to obtain the frame bias components, the routine `P06E` to obtain various precession
angles, and the routine `NUM06A` to generate the nutation matrix. The product N×P×B is formed
by using the routine `RXR` ([^2] this is matrix product). Next call the routine `GST06` to obtain the GST corresponding to the
current UT1, and apply it as an R3 rotation using the routine RZ to form the true matrix-to-TIRS.
Given xp, yp, and obtaining s_0 with the routine `SP00`, the polar motion matrix (i.e. TIRS-to-ITRS)
is then obtained using the routine `POM00`. The product of the two matrices (GCRS-to-TIRS and
TIRS-to-ITRS), obtained by calling the routine `RXR` ([^2]), is the GCRS-to-ITRS matrix, which can be
inverted by calling the routine `TR` ([^2]) to give the final result.

Methods (1) and (2) agree to microarcsecond precision.

The abridged nutation model IAU 2000B (see Section 5.5.1) can be substituted in Method (2) by
calling `NUT00B` instead of `NUT06A`. Depending on the application, the best compromise between
speed and accuracy may be to evaluate the full series to obtain sample values for interpolation.
