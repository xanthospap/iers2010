# libiers10++

[![Build Status](https://app.travis-ci.com/xanthospap/iers2010.svg?branch=master)](https://app.travis-ci.com/xanthospap/iers2010)

C++ library implementing the IERS 2010 standards.

# Introduction

This project contains a number of functions implementing models defined in
[IERS Conventions (2010)](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html/) 
as described in [^2].
The International Earth Rotation and Reference Systems Service ([IERS](https://www.iers.org/IERS/EN/Home/home_node.html))
publishes the Conventions along with relevant documents, model implementations and 
respective test cases; the latter two are available in the FORTRAN programming 
language at the [IERS website](https://iers-conventions.obspm.fr/conventions_material.php).

This repository is an effort to translate the algorithms in the C++ programming 
language with (as much as possible) minor modifications. Note that the software 
found at this website is routinely updated.

# Prerequisites <a name="prerequisites"></a>

The C++ library [ggdatetime](https://github.com/xanthospap/ggdatetime) 
is used in the library to handle datetime instances when needed. 
Hence, you should have [ggdatetime](https://github.com/xanthospap/ggdatetime) 
on your system.

The C++ library [ggeodesy](https://github.com/xanthospap/ggeodesy) is used 
to handle basic geodesy when needed. Again, you should have the library 
installed (on your system) to successefuly compile this library.

Vector/Matrix operations are implemented using the [eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page) 
library. Availability of this library is also mandatory.

# CSPICE <a name="cspice"></a>

Interaction with planetary ephemeris files (i.e. JPL's published DE), is 
done via the [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html). 
The toolkit is a C library with a corresponding API. 

## Installation

* Download the C library from the [correponding site](https://naif.jpl.nasa.gov/naif/toolkit_C.html) and uncompress.

* Use the script [c2cpp_header.py](script/cppspice/c2cpp_header.py) to tranform 
C header file. Run the script using the cspice `include` folder path as command 
line argument. I.e. if the uncompressed cspice folder is at `/home/work/var/cspice`, 
use `$> c2cpp_header.py /home/work/var/cspice`.

* Run the `makeall.csh` script provided by the distribution (under `cspice` folder). Note 
that the script is in the C-sheel, hence you might need to `$>csh makeall.csh`.

* Copy the [install.sh](script/cppspice/install.sh) script under the `cspicer` folder; 
run it as root, to install the library. Header files will be available at 
`/usr/local/include` and the library at `/usr/local/lib`.

# Compilation / Installation <a name="compilation-installation"></a>

> January 2022:
> From now on, only the [scons](https://scons.org/) build will be supported; 
> support for autotools is dropped.

The project is built via the [scons](https://scons.org/) built system:

```bash
$> git clone https://github.com/xanthospap/iers2010.git && cd iers2010
$> scons
$> sudo scons install
```

## Test Programs <a name="test-programs"></a>
To compile the test utilities, use the `make-test` switch (i.e. `scons make-test=1`). 
Note that the SOFA [^1] library is needed to compile and run the tests.

# Library Components <a name="library-components"></a>

## Earth Orientation Parameters (EOPs)

The library can parse EOP information from IERS files using the IERS C04/14 
(see [updateC04.txt](https://hpiers.obspm.fr/iers/eop/eopc04_14/updateC04.txt)) and 
IERS C04/20 (see [eopc04.txt](https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.txt)) 
file formats/series. The files can be downloaded from IERS, e.g. 
[eopc04_IAU2000.62-now](https://hpiers.obspm.fr/iers/eop/eopc04_14/eopc04_IAU2000.62-now) 
and [eopc04.1962-now](https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now).

Normally, EOP information is stored in an `EopSeries` instance, which can hold 
records for multiple epochs (ordered chronologically) and perform basic utilities 
such as interpolation.

Users can choose to "regularize" ΔUT1 and LOD values, i.e. remove zonal tidal 
variations with frequencies ranging from 5 days to 8.6 years from UT1 values 
(result typically denoted UT1R) and LOD (see Chapter 8.1 of [^2] and [^3]).

The most important and usual operation an `EopSeries` instance performs, is the 
interpolation of EOP values for an epoch of interest. This can be done with a 
simple call to `EopSeries::interpolate` method. Note that for utmost accuracy 
certain corrections will have to be applied to the resulting interpolated values, 
such as removal of libration and ocean tidal effects (see the example source code).

For a complete example of `EopSeries` usage, see [eop_interpolation.cpp](examples/eop_interpolation.cpp).

## Notes on Dates and Datetimes <a name="notes-on-dates"></a>
This library makes use of [ggdatetime](https://github.com/xanthospap/ggdatetime) 
for handling dates and datetime instances. By default, this means that we are 
(mostly) handling datetimes in the so-called **Modified Julian Date** format, 
instead of the **Julian Date** format used throughout most of the SOFA [^1] functions. 

## Astronomy Functionality and Celestial Frames <a name="astronomy"></a>
The library contains a list of functions for basic astronomical computations, 
especially related to computations involved in transforming from/to the 
Geocentric Celestial Reference Frame (GCRF). Declerations can be found in the 
header files [iau.hpp](src/iau.hpp) and [fundarg.hpp](src/fundarg.hpp).

There are a number of ways to transform between (geocentric) inertial and earth-fixed 
frames; in this library, we target the so called **"CIO-based"** transformation, based 
on the **IAU 200/2006** model (see [^1] and [^2]), as this is the recommended approach 
as described in the IERS2010 standards.

## Tests and Compliance with IAU/SOFA <a name="tests-and-compliance"></a>

Once you have [compiled the test-suite](#test-programs), you can use the script 
`run_sofa_checks.py` to check the compliance with respect to the IAU SOFA 
library (e.g. `run_sofa_checks.py --progs-dir test/sofa_unit_tests/`). This 
will run a number of checks and provide discrepancies (w.r.t. SOFA) for all 
functions/parameters common to both libraries. Explicit results can be found in
[this table](#sofa-check-table).

### Checking Descripancies <a name="checking-descripancies"></a>

Descripancies for angular parameters are usually checked in units of *arcseconds*.
Descripancies for rotation matrices (.e.g. $R_A$ and $R_B$) are obtained by the 
*axis-angle* of the combined matrix as:
$\theta = arccos(\frac{(R_A^T R_B) -1}{2})$ again transformed in *arcseconds*.

### Result Key-Points

* Discrepancies between the SOFA implementation and `libiers2010` for one day of 
GRACE orbit are presented in ??. They ammount up to $\approx 1\times10^{-4} mm$ 
for the $X$ and $Y$ components and $\approx 2\times10^{-6} mm$ for $Z$. See results 
[here](#sofa_ter2cel-dX)

* The accuracy of the transformation is presented in presented in ??. It is 
$\approx 5\times10^{-6} mm$ for each of the $X$, $Y$ and $Z$ components.(Note 
that this is checked again using one day of GRACE orbit, and performing the 
circular transformation `ICRF->ITRF->ICRF`). See results [here](#internal_ter2cel-dX)

## Tests and Compliance with IERS Implementation <a name="tests-and-compliance-fortran"></a>

Once you have [compiled the test-suite](#test-programs), you can use the script 
`run_fortran_checks.py` to check the compliance of the current library against 
the FORTRAN subroutines distributed by the IERS.
Example run:
```
$>./run_fortran_checks.py --for-dir fortran_impl/ \
    --cpp-dir test/fortran_unit_tests/ \
    --verbose
```
To run these unit tests however, you will need to have downloaded the IERS-distributed 
source code accompanying the standards. Explicit results can be found in 
[this table](#fortran-check-table).

### Computation of GMST

In this implementation computation of GMST is performed using the IAU2006 
model, i.e. using both UT1 and TT (see Sec. 5.5.7, ERA based expressions for 
Greenwich Sidereal Time). Contrary to this, some IERS-distributed FORTRAN 
subroutines use an older model for GMST (i.e. `FCNNUT`, `PMSDNUT2`).

# References <a name="references"></a>

[^1]: International Astronomical Union, Standards of Fundamental Astronomy (SOFA) 
software tools collection, available at [https://www.iausofa.org/]

[^2]: Gérard Petit and Brian J. Luzum, editors. IERS Conventions (2010), volume 36 of
IERS Technical Note, 2010. International Earth Rotation and Reference Systems
Service (IERS), International Earth Rotation and Reference Systems Service
(IERS).

[^3] Bradley, B. K., A. Sibois, and P. Axelrad (2016). “Influence of ITRS/GCRS im-
plementation for astrodynamics: Coordinate transformations”. In: Advances in
Space Research 57.3, pp. 850–866. issn: 0273-1177. doi: https://doi.org/
10.1016/j.asr.2015.11.006. url: https://www.sciencedirect.com/
science/article/pii/S0273117715007929.

# Acknowledgement

Software Routines from the IAU SOFA Collection were used. Copyright © International Astronomical Union Standards of Fundamental Astronomy (http://www.iausofa.org)
