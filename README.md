# libiers10++

[![Build Status](https://app.travis-ci.com/xanthospap/iers2010.svg?branch=master)](https://app.travis-ci.com/xanthospap/iers2010)

C++ library implementing the IERS 2010 standards.

Table of Contents
1. [Introduction](#introduction)
2. [Prerequisites](#introduction)
3. [Compilation / Installation](#installation)
4. [Data Files](#data-files)
5. [BLQ format files](#blq-files)
6. [Status of the Library _(or translated routines)_](#library-status)
7. [Available functions/models **not** part of the IERS 2010 Conventions](#library-extra)
8. [Test Cases & FORTRAN vs C++ implementations](#test-cases)
1. [Test Programs _(optional read)_](#test-programs)
9. [How to use the library](#how-to)
10. [Documentation & Library API (TODO)](#dox)
11. [TODO](#todo)
12. [Bugs & Maintanance](#bugs)

# Introduction <a name="introduction"></a>
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
for the $X$ and $Y$ components and $\approx 2\times10^{-6} mm$ for $Z$.

* The accuracy of the transformation is presented in presented in ??. It is 
$\approx 5\times10^{-6} mm$ for each of the $X$, $Y$ and $Z$ components.(Note 
that this is checked again using one day of GRACE orbit, and performing the 
circular transformation `ICRF->ITRF->ICRF`)

# References <a name="references"></a>

[^1]: International Astronomical Union, Standards of Fundamental Astronomy (SOFA) 
software tools collection, available at [https://www.iausofa.org/]

[^2]: Gérard Petit and Brian J. Luzum, editors. IERS Conventions (2010), volume 36 of
IERS Technical Note, 2010. International Earth Rotation and Reference Systems
Service (IERS), International Earth Rotation and Reference Systems Service
(IERS).

# Acknowledgement

Software Routines from the IAU SOFA Collection were used. Copyright © International Astronomical Union Standards of Fundamental Astronomy (http://www.iausofa.org)

# Results <a name="results"></a>

## Discrepancies w.r.t IAU/SOFA <a name="sofa-check-table"></a>

function  |argument  |num tests |num fails |max error |param. type |status    
----------|----------|----------|----------|----------|------------|----------
nut00a    |dpsi      |    100000|     49379|-1.202e-12|Angle       |FAILED    
nut00a    |deps      |    100000|     48827|4.920e-13|Angle       |FAILED    
pom00     |          |    100000|         0|0.000e+00|RotMatrix   |OK        
s00       |          |    100000|         0|0.000e+00|Angle       |OK        
pn00a     |dpsi      |    100000|     49377|-1.201e-12|Angle       |FAILED    
pn00a     |deps      |    100000|     48545|4.827e-13|Angle       |FAILED    
pn00a     |epsa      |    100000|         0|0.000e+00|Angle       |OK        
pn00a     |rb        |    100000|         0|0.000e+00|RotMatrix   |OK        
pn00a     |rp        |    100000|     17689|7.529e-03|RotMatrix   |FAILED    
pn00a     |rbp       |    100000|     20097|7.529e-03|RotMatrix   |FAILED    
pn00a     |rn        |    100000|     53957|6.147e-03|RotMatrix   |FAILED    
pn00a     |rbpn      |    100000|     36207|9.720e-03|RotMatrix   |FAILED    
num06a    |          |    100000|     42829|6.147e-03|RotMatrix   |FAILED    
bi00      |dpsibi    |    100000|         0|0.000e+00|Angle       |OK        
bi00      |depsbi    |    100000|         0|0.000e+00|Angle       |OK        
bi00      |dra       |    100000|         0|0.000e+00|Angle       |OK        
c2ixys    |          |   1000000|         0|0.000e+00|RotMatrix   |OK        
pnm06a    |          |    100000|     44542|7.529e-03|RotMatrix   |FAILED    
bp00      |rb        |    100000|         0|0.000e+00|RotMatrix   |OK        
bp00      |rp        |    100000|     17677|7.529e-03|RotMatrix   |FAILED    
bp00      |rbp       |    100000|     20188|7.529e-03|RotMatrix   |FAILED    
fal03     |fapa03    |    100000|         0|0.000e+00|Angle       |OK        
fal03     |fane03    |    100000|         0|0.000e+00|Angle       |OK        
fal03     |faur03    |    100000|         0|0.000e+00|Angle       |OK        
fal03     |fasa03    |    100000|         0|0.000e+00|Angle       |OK        
fal03     |faju03    |    100000|         0|0.000e+00|Angle       |OK        
fal03     |fama03    |    100000|         0|0.000e+00|Angle       |OK        
fal03     |fae03     |    100000|         0|0.000e+00|Angle       |OK        
fal03     |fave03    |    100000|         0|0.000e+00|Angle       |OK        
fal03     |fame03    |    100000|         0|0.000e+00|Angle       |OK        
fal03     |faom03    |    100000|         0|0.000e+00|Angle       |OK        
fal03     |fad03     |    100000|         0|0.000e+00|Angle       |OK        
fal03     |faf03     |    100000|         0|0.000e+00|Angle       |OK        
fal03     |falp03    |    100000|         0|0.000e+00|Angle       |OK        
fal03     |fal03     |    100000|         0|0.000e+00|Angle       |OK        
fw2m      |          |    100000|         0|0.000e+00|Angle       |OK        
sp00      |ee06a     |    100000|     62460|1.832e-10|Angle       |FAILED    
sp00      |gmst06    |    100000|      4612|1.172e-08|Angle       |FAILED    
sp00      |gmst00    |    100000|      4603|1.172e-08|Angle       |FAILED    
sp00      |gst06a    |    100000|      4678|5.863e-09|Angle       |FAILED    
sp00      |obl06     |    100000|         2|2.290e-11|Angle       |FAILED    
sp00      |obl80     |    100000|         0|0.000e+00|Angle       |OK        
sp00      |ee00      |    100000|         0|0.000e+00|Angle       |OK        
sp00      |eect00    |    100000|     25677|2.150e-17|Angle       |FAILED    
sp00      |sp00      |    100000|      8006|5.332e-21|Angle       |FAILED    
era00     |          |    100000|      2404|1.172e-08|Angle       |FAILED    
xy06      |xp        |     10000|      1528|-6.262e-13|Angle       |FAILED    
xy06      |yp        |     10000|      3272|-4.626e-13|Angle       |FAILED    
nut06a    |dpsi      |     10000|      4926|-1.037e-12|Angle       |FAILED    
nut06a    |deps      |     10000|      4824|-4.548e-13|Angle       |FAILED    
c2i06a    |          |     10000|      5392|6.147e-03|RotMatrix   |FAILED    
pn00      |rb        |     10000|         0|0.000e+00|RotMatrix   |OK        
pn00      |rp        |     10000|      1821|6.147e-03|RotMatrix   |FAILED    
pn00      |rbp       |     10000|      2060|7.529e-03|RotMatrix   |FAILED    
pn00      |rn        |     10000|         0|0.000e+00|RotMatrix   |OK        
pn00      |rbpn      |     10000|       456|7.529e-03|RotMatrix   |FAILED    
pn00      |epsa      |     10000|         0|0.000e+00|Angle       |OK        
numat     |          |     10000|         0|0.000e+00|RotMatrix   |OK        
pr00      |dpsipr    |     10000|       621|2.184e-17|Angle       |FAILED    
pr00      |depspr    |     10000|       108|2.730e-18|Angle       |FAILED    
pfw06     |gamb      |     10000|       978|-1.398e-15|Angle       |FAILED    
pfw06     |phib      |     10000|         0|0.000e+00|Angle       |OK        
pfw06     |psib      |     10000|       878|-5.367e-13|Angle       |FAILED    
pfw06     |epsa      |     10000|         0|0.000e+00|Angle       |OK        
pn06      |rb        |     10000|     10000|6.147e-03|RotMatrix   |FAILED    
pn06      |rp        |     10000|     10000|8.693e-03|RotMatrix   |FAILED    
pn06      |rbp       |     10000|      3952|7.529e-03|RotMatrix   |FAILED    
pn06      |rn        |     10000|      1217|8.693e-03|RotMatrix   |FAILED    
pn06      |rbpn      |     10000|      1369|6.147e-03|RotMatrix   |FAILED    
pn06      |epsa      |     10000|         0|0.000e+00|Angle       |OK        
p06e      |eps0      |    100000|         0|0.000e+00|Angle       |OK        
p06e      |psia      |    100000|      8446|7.156e-13|Angle       |FAILED    
p06e      |oma       |    100000|         0|0.000e+00|Angle       |OK        
p06e      |bpa       |    100000|      7004|6.989e-16|Angle       |FAILED    
p06e      |bqa       |    100000|      8284|-5.591e-15|Angle       |FAILED    
p06e      |pia       |    100000|      8149|-5.591e-15|Angle       |FAILED    
p06e      |bpia      |    100000|         4|1.832e-10|Angle       |FAILED    
p06e      |epsa      |    100000|         0|0.000e+00|Angle       |OK        
p06e      |chia      |    100000|      8791|-1.398e-15|Angle       |FAILED    
p06e      |za        |    100000|      7996|3.578e-13|Angle       |FAILED    
p06e      |zetaa     |    100000|      8049|-3.578e-13|Angle       |FAILED    
p06e      |thetaa    |    100000|      6403|1.789e-13|Angle       |FAILED    
p06e      |pa        |    100000|      8500|-7.156e-13|Angle       |FAILED    
p06e      |gam       |    100000|      8820|1.398e-15|Angle       |FAILED    
p06e      |phi       |    100000|         0|0.000e+00|Angle       |OK        
p06e      |psi       |    100000|      8636|-7.156e-13|Angle       |FAILED    
c2t06a    |          |     10000|      3042|7.529e-03|RotMatrix   |FAILED    
xys06a    |xp        |     10000|      1986|-5.367e-13|Angle       |FAILED    
xys06a    |yp        |     10000|      2402|-3.435e-11|Angle       |FAILED    
xys06a    |s         |     10000|      4329|4.703e-14|Angle       |FAILED    
xys00a    |xp        |     10000|      2381|6.262e-13|Angle       |FAILED    
xys00a    |yp        |     10000|       417|-1.145e-11|Angle       |FAILED    
xys00a    |s         |     10000|      3324|1.687e-14|Angle       |FAILED    
s06       |          |    100000|         0|0.000e+00|Angle       |OK
<!--
# Data Files <a name="data-files"></a>

A couple of (data) files are needed to use some functions wihtin the library, 
namely:

  * the data file `gpt2_5.grd` is needed for the function `gpt2` (distributed by 
  [VMF Data Server](https://vmf.geo.tuwien.ac.at/) at TU Wien [^1]), available 
  at https://vmf.geo.tuwien.ac.at/codes/gpt2_5.grd

  * either (or both) of the files `gpt3_5.grd` and `gpt3_1.grd` () distributed by 
  [VMF Data Server](https://vmf.geo.tuwien.ac.at/) at TU Wien [^1]), available 
  at https://vmf.geo.tuwien.ac.at/codes/gpt3_5.grd and 
  https://vmf.geo.tuwien.ac.at/codes/gpt3_1.grd . They are needed to use functions 
  related to the GPT3 tropospheric mapping functions.

[^1]: VMF Data Server; editing status 2020-12-14; re3data.org - Registry of Research 
Data Repositories. http://doi.org/10.17616/R3RD2H

~~This repository includes also the data file [gpt2_5.grd](data/gpt2_5.grd) in the
`data` directory. This file is needed for computations when the function [gpt2](#gpt2-cmp) in invoked.
When installing the libary (aka `make install`) this file will be installed in the 
default `share` directory of the system, under the folder `iers10`. In most *X systems, 
this means that you'll end up with the file: `/usr/local/share/iers10/gpt2_5.grd`.~~

~~If you want to change the data directory path, you will need to alter the respective 
Makefile template, that is [data/Makefile.am](data/Makefile.am).~~

# BLQ format files <a name="blq-files"></a>

The library includes a helper class, i.e. `iers2010::BlqIn` (declared in 
[blqstrm.hpp](src/blqstrm.hpp)) to assist reading records off from a BLQ file. 
Obviously, users can make use of this code independent of the (rest of the) 
library.

The file [test_blq.cpp](test/test_blq.cpp) includes a test case usage of the 
source code for reading and manipulating BLQ files; the source code is compiled to 
the executable `testBlq`.

# Status of the Library _(or translated routines)_ <a name="library-status"></a>

| Chapter | (Sub)Routine | Translated | Tested | Version  | Comments |
|:--------|:-------------|:----------:|:------:|:---------|:---------|
| 4       | [GCONV2](https://iers-conventions.obspm.fr/content/chapter4/GCONV2.F)              |<ul><li>- [ ] </li></ul>|<ul><li>- [ ] </li></ul>| | see [ngpt car2ell](https://github.com/xanthospap/ngpt/blob/master/src/car2ell.hpp)|
| 5       | [PMSDNUT2](https://iers-conventions.obspm.fr/content/convupdt/chapter5/PMSDNUT2.F) |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 13.11.2011 | see [pmsdnut2](#pmsdnut2-cmp) Discrepancies between FORTRAN and C++ implementation <1e<sup>-6</sup> mas; this is due to the D0 decleration|
|         | [UTLIBR](https://iers-conventions.obspm.fr/content/chapter5/UTLIBR.F)              |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 23.06.2010 | see [utlibr](#utlibr-cmp) Discrepancies between FORTRAN and C++ implementation ~1e-7 or less; this is due to the     D0 decleration|
|         | [FUNDARG](https://iers-conventions.obspm.fr/content/chapter5/FUNDARG.F)            |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 25.02.2010 | see [fundarg](#fundarg-cmp) C++ and FORTRAN implementations give identical results,which do not agree though with the test case provided |
|         | [FCNNUT](https://iers-conventions.obspm.fr/content/convupdt/chapter5/FCNNUT.F)     |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 19.12.2013 | |
| 7       | [DEHANTTIDEINEL](https://iers-conventions.obspm.fr/content/convupdt/chapter7/dehanttideinel/DEHANTTIDEINEL.F) |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 19.12.2016 | see [dehanttideinel](#dehanttideinel-cmp). One test case does not produce the expected results both in the FORTRAN and the C++ implementation; that is test case 4, as privided in the source code|
|         | [HARDISP](https://iers-conventions.obspm.fr/content/convupdt/chapter7/hardisp/HARDISP.F) |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 19.12.2016 | see [hardisp](#hardisp-cmp) |
|         | [ARG2](https://iers-conventions.obspm.fr/content/convupdt/chapter7/ARG2.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 07.10.2011 | see [arg2](#arg2-cmp)|
|         | [IERS_CMP_2015](https://iers-conventions.obspm.fr/content/convupdt/chapter7/IERS_CMP_2015.F) | <ul><li>- [ ] </li></ul>|<ul><li>- [ ] </li></ul>| | |
| 8       | [CNMTX](https://iers-conventions.obspm.fr/content/chapter8/CNMTX.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 17.03.2010 | see [cnmtx](#cnmtx-cmp) |
|         | [ORTHO_EOP](https://iers-conventions.obspm.fr/content/chapter8/ORTHO_EOP.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 19.05.2010 | see [ortho_eop](#ortho_eop-cmp) |
|         | [RG_ZONT2](https://iers-conventions.obspm.fr/content/convupdt/chapter8/RG_ZONT2.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 20.12.2011 | see [rg_zont2](#rg_zont2-cmp) |
| 9       | [FCUL_A](https://iers-conventions.obspm.fr/content/chapter9/FCUL_A.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 13.08.2009 | see [fcul_a](#fcul_a-cmp) |
|         | [FCUL_B](https://iers-conventions.obspm.fr/content/chapter9/FCUL_B.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 14.08.2009 | see [fcul_b](#fcul_b-cmp) |
|         | [FCUL_ZD_HPA](https://iers-conventions.obspm.fr/content/chapter9/FCUL_ZD_HPA.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 14.08.2009 | see [fcul_zd_hpa](#fcul_zd_hpa-cmp) Both FORTRAN and C++ implementations prodice different results from the test case provided|
|         | [FCUL_ZD_HPA](https://iers-conventions.obspm.fr/content/chapter9/FCUL_ZD_HPA.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 14.08.2009 | see [fculzd_hpa](#fculzd_hpa-cmp) |
|         | [GMF](https://iers-conventions.obspm.fr/content/chapter9/GMF.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 12.08.2009 | see [gmf](#gmf-cmp) |
|         | [VMF1](https://iers-conventions.obspm.fr/content/chapter9/VMF1.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 12.01.20012 | see [vmf1](#vmf1-cmp) |
|         | [VMF1_HT](https://iers-conventions.obspm.fr/content/chapter9/VMF1_HT.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 12.01.20012 | see [vmf1_ht](#vmf1_ht-cmp) |
|         | [GPT](https://iers-conventions.obspm.fr/content/chapter9/GPT.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 18.10.2011 | see [gpt](#gpt-cmp) |
|         | [GPT2](https://iers-conventions.obspm.fr/content/chapter9/GPT2.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 31.05.2013 | see [gpt2](#gpt2-cmp) |


# Available functions/models **not** part of the IERS 2010 Conventions <a name="library-extra"></a>

The library includes a number of functions and/or model implementations that are 
not part of the the IERS 2010 Conventions, but are widely used in the field of 
Satellite Geodesy. They include:

  * the IERS published function [INTERP](https://hpiers.obspm.fr/iers/models/interp.f)
  used to ocean tidal and libration effects in pole coordinates
  
  * a number functions implementing models of Earth motion and orientation, 
  astrometry and Celestial to Terrestrial RF transformations. These follow the IAU 
  Resolutions (2000/2006) and are included via (the header file) [iau.hpp](src/iau.hpp).
  A list and details are recorded in [readme.md](src/iau/readme.md).
  
  * a number functions implementing atmospheric models models, listed in the 
  directory [extra/atmosphere](src/extra/atmosphere/)

# Test Cases & FORTRAN vs C++ implementations <a name="test-cases"></a>

> The discussion below is only a description of the discrepancies between implementations
> and does **not** depict the actual precission of the models.

### Test Programs _(optional read)_ <a name="test-programs"></a>

To build the tests, enter: `scons --make-check`. The test programs will be build 
in the `test` directory, with names `test-[function name].out`.

| Function | Test Case | Comments |
|:-------- |:----------|:---------|
| fcnnut   | Use [test_fcnnut.f](fortran_impl/test_fcnnut.f) to compute results using the original FORTRAN implementation. Copy to [test_fcnnut.cpp](test/test_fcnnut.cpp)| Dispcrepancies < 1e-9 [μas] |
| fundarg  | Use [test_fundarg.f](fortran_impl/test_fundarg.f) to compute results using the original FORTRAN implementation. Copy to [test_fundarg.cpp](test/test_fundarg.cpp)| Dispcrepancies < 1e-11 [rad] |
| pmsdnut2 | Use [test_pmsdnut2.f](fortran_impl/test_pmsdnut2.f) to compute results using the original FORTRAN implementation. Copy to [test_pmsdnut2.cpp](test/test_pmsdnut2.cpp)| Dispcrepancies < 1e-6 [μas] [^pmsdnut2]|
| utlibr | Use [test_utlibr.f](fortran_impl/test_utlibr.f) to compute results using the original FORTRAN implementation. Copy to [test_utlibr.cpp](test/test_utlibr.cpp)| Dispcrepancies [^utlibr] <ul><li> < 1e-7 [μas] for dut1 and</li><li> < 1e-5 [μas/day] for dlod</li></ul> |
| arg2  | Use [test_arg2.f](fortran_impl/test_arg2.f) to compute results using the original FORTRAN implementation. Copy to [test_arg2.cpp](test/test_arg2.cpp)| Dispcrepancies < 1e-10 radians |
| oeop::cnmtx | Use [make_test_cnmtx](fortran_impl/make_test_cnmtx.f) and copy output to [test_cnmtx_random.cpp](test/test_cnmtx_random.cpp) | Note that the function is used within `ortho_eop`. Results agree to better than 1e-12 [-] |
| ortho_eop  | Use [test_orthoeop.f](fortran_impl/test_orthoeop.f) to compute results using the original FORTRAN implementation. Copy to [test_orthoeop.cpp](test/test_orthoeop.cpp)| Dispcrepancies < 1e-12 [μas] |
| dehanttideinel | Use [dehanttideinel/MAIN.F](fortran_impl/dehanttideinel/MAIN.F) to compute results using the original FORTRAN implementation. Copy to [test_dehanttideinel.cpp](test/test_dehanttideinel.cpp); see [readme](fortran_impl/dehanttideinel/readme.md) in the folder | Dispcrepancies < 1e-12 [m]. Note that this routine is altered quite heavily, towards computing displacements for a number of sites, at the same instant in time. |

[^pmsdnut2]: Dispcrepancies < 1e-9 μas can be obtained if we change the declerations 
in the original FORTRAN source code, in the `DATA` matrix. Double-precision numerics, 
can be extended with `D0` (e.g. `6798.3837D0,  0.0D0,  0.6D0, -0.1D0, -0.1D0,`, 
instead of `6798.3837, 0.0,   0.6,   -0.1,   -0.1,`)

[^utlibr]: Dispcrepancies < 1e-9 [μas] and 1e-8 [μas/day] can be obtained if we change the declerations 
in the original FORTRAN source code, in the `DATA` matrix.

~~To compile the test programs, you need to enter the command `make check` at the 
`ROOTDIR` folder (after you have run `make`). This will build the programs 
to test the implementations of the individual functions in the library. 
They are compiled into executables in the `test` folder:~~
  
  * ~~[testHardisp](test/test_hardisp.cpp) compiled to `test/testHardisp`
  This program checks the library hardisp program; to run this, you are going to need
  a `BLQ` file with records for the GNSS station 'ONSA' (or you could use the file 
  [NTUA.BLQ](data/NTUA.BLQ)). For details see [hardisp](#hardisp-cmp).~~

- ~~[testIers2010](test/test_iers2010.cpp) compiled to `test/testIers2010`
  This program checks the library functions using the test cases provided in the original
  FORTRAN routines. Run and check the results; in some cases, the differences (if 
  present) are expected and reason they appear is documented in the routine-specific
  blocks below.~~

- ~~[testDehantTide](test/test_dehanttide.cpp) compiled to `test/testDehantTide`
  This program checks specifically the library function `dehandtideinel` using the test cases 
  provided in the original FORTRAN routine.~~

~~These programs check the library functions against the test cases provided in 
the IERS published FORTRAN source code files. **Note that in some rare occasions, the 
original FORTRAN implementation may fail (to a given accuracy) the actual test case** 
(see e.g. [fundarg](#fundarg-cmp))~~

~~If needed, alternative FORTRAN implementations are provided (in the `fortran_impl` directory) for further testing.
These can be compiled using the `fortran_impl/Makefile` aka run `make` in the `fortran_impl` folder.
Note that you will need the IERS source code (you can use the script [fetch_iers_lib.sh](fortran_impl/fetch_iers_lib.sh))
and [gfortran](https://gcc.gnu.org/fortran/). Two programs are compiled, `fortran_impl/test_iers` and
`fortran_impl/test_iers_d0`; the latter uses the alternative FORTRAN implementations (**NOT** the original
IERS source code) Alternative FORTRAN implementations are provided for:~~

~~-PMSDNUT2.F named [PMSDNUT2_D0.F](fortran_impl/PMSDNUT2_D0.F) (see [pmsdnut2](#pmsdnut2-cmp))~~
~~-UTLIBR.F named [UTLIBR_D0.F](fortran_impl/UTLIBR_D0.F) (see [utlibr](#utlibr-cmp))~~

### fundarg <a id="fundarg-cmp"></a>

Test case provided in the `FUNDARG.F` source file; the C++ implementation 
([test_fundarg.cpp](src/test_fundarg.cpp)) and the FORTRAN implementation 
([test_fundarg.f](fortran_impl/test_fundarg.f)) produce the same descripancies 
when compared to the test case given; these are:
```
----------------------------------------
> fundarg
----------------------------------------

args[0] = 7.638334e-14 radians
args[1] = 9.699797e-12 radians
args[2] = 7.771561e-14 radians
args[3] = 9.591439e-12 radians
args[4] = 2.220446e-16 radians
```

### pmsdnut2 <a id="pmsdnut2-cmp"></a>

Test case provided in the `PMSDNUT2.F` source file; the C++ implementation shows 
discrepancies of
```
----------------------------------------
> pmsdnut2
----------------------------------------

dx= 1.328753e-07 microarcseconds
dy= 4.441797e-07 microarcseconds
```
These discrepancies are due to the initialization of double precision values in 
the PMSDNUT2.F subroutine; if we initialize these values using the notation:
`123.D0` instead of `123.0`, results are identical for the FORTRAN and C++ 
implementations.

*Note that the Same thing happens with [utlibr](#utlibr-cmp).*

### utlibr <a id="utlibr-cmp"></a>

The source file UTLIBR.F includes two individual test cases; the C++ 
implementation shows discrepancies of
```
----------------------------------------
> utlibr
----------------------------------------

dut1= 1.373660e-08 mas
dlod= 2.357558e-07 mas / day
dut1= 1.522151e-08 mas
dlod= 4.876766e-07 mas / day
```

This is because of the initialization of the DOUBLE PRECISION float arrays 
`PER, DUT1S, DUT1C, DLODS and DLODC` which not explicitely marked with `D0` 
(in `UTLIBR.F`). If we change the initialization, the C++ and FORTRAN versions 
produce identical results.

*Same thing happens with [pmsdnut2](#pmsdnut2-cmp).*

### fcnnut <a id="fcnnut-cmp"></a>

FORTRAN (i.e. `FCNNUT.F`) and C++ implementations produce identical results.

### arg2 <a id="arg2-cmp"></a>

FORTRAN (i.e. `ARG2.F`) and C++ implementations produce identical results.

### dehanttideinel <a id="dehanttideinel-cmp"></a>

There are four test cases available (in `DEHANTTIDEINEL.F`) for this routine. 
For the first three, they all yield discrepancies < 1e<sup>-13</sup> meters for the C++ 
implementation, while for the FORTRAN implementation differences are zero. __However, the 
fourth test case fails misserably, both for the C++ and the FORTRAN implementation__; in this 
case, discrepancies in the order of tens of meters are found! I do not know what 
the problem is.

### cnmtx <a id="cnmtx-cmp"></a>

FORTRAN (i.e. `CNMTX.F`) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### ortho_eop <a id="ortho_eop-cmp"></a>

FORTRAN (i.e. `ORTHO_EOP.F`) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### rg_zont2 <a id="rg_zont2-cmp"></a>

The test case (provided in `RG_ZONT2.F`) implemented using the C++ function, gives max discrepancies 
in the order of 1e<sup>-19</sup>.

### fcul_a <a id="fcul_a-cmp"></a>

FORTRAN (i.e. `FCUL_A.F`) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### fcul_b <a id="fcul_b-cmp"></a>

FORTRAN (i.e. `FCUL_B.F`) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### fcul_zd_hpa <a id="fcul_zd_hpa-cmp"></a>

FOFTRAN and C++ implementations produce identical results; however, these show 
discrepancies up to 5e<sup>-6</sup> with the ones provided at the test case (in the header of `FZUL_ZD_HPA.F`).

### gmf <a id="gmf-cmp"></a>

FORTRAN (i.e. `GMF.F`) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### vmf1 <a id="vmf1-cmp"></a>

FORTRAN (i.e. `VMF1.F`) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### vmf1_ht <a id="vmf1_ht-cmp"></a>

FORTRAN (i.e. `VMF1_HT.F`) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### gpt <a id="gpt-cmp"></a>

FORTRAN (i.e. `GPT.F`) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### gpt2 <a id="gpt2-cmp"></a>

FORTRAN (i.e. `GPT2.F`) and C++ implementations produce identical results. The test case 
provided in the FORTRAN source code, is limited to the accuracy actually provided by the 
algorithm/model; to that accuracy, both implementations produce identical results.

### hardisp <a id="hardisp-cmp"></a>

There is currently a problem with HARDISP algorithm; a date in UTC is given as input which is then
converted to ET (aka Ephemeris Time) and used to compute results. However, according to IAU, ET
is not a time-scale that should be used and is actually not supported (see International Astronomical Union, 
Standards Of Fundamental Astronomy, SOFA Time Scale and Calendar Tools, Software version 13 - 
Document revision 1.6). In this document, it is stated that ``ET (ephemeris time): superseded by TT and TDB``

In the C++ implementation instead of ET we use TT. To do that, we use the formula: 
``UTC + ΔAT + 32.184(sec) = TT``. The transformation (along with other date/time
computations) are performed in the function tdfrph. The two implementations (aka
the FORTRAN and C++) produce differences up to 0.000001 meters.

To test hardisp (the standalone program that is), use e.g.
```
src/hardisp 2009 6 25 1 10 45 24 3600 < test/onsa.blq > onsa.tmp
paste onsa.tmp test/test_onsa_results | awk '{printf "%9.6f %9.6f %9.6f\n", $1-$4, $2-$5, $3-$6}'
```
OR
```
src/hardisp 2009 6 25 1 10 45 24 3600 < test/reyk.blq > reyk.tmp
paste reyk.tmp test/test_reyk_results | awk '{printf "%9.6f %9.6f %9.6f\n", $1-$4, $2-$5, $3-$6}'
```

You should see the differences (per line and column) and they should not exceed
the value 1e<sup>-6</sup>.

In the official IERS2010 software, HARDISP comes as a standalone routine/program. In
this implementation, HARDISP comes both as a standalone program as well as a function 
(see `hisp::hardisp_impl` in file [hardisp_impl.cpp](src/hardisp_impl.cpp)) that can be 
used by the users in source code.

A test case for using the implementation is provided in the test source code 
[test_hardisp.cpp](test/test_hardisp.cpp). When compiling the tests (via `make check`) 
this source code will be compiled to `test/testHardisp`. Given a BLQ file with a 
record for station ONSA, you can check the results (which sould exactly agree with 
the `src/hardisp` executable. A BLQ file with records for ONSA is provided in [NTUA.BLQ](data/NTUA.BLQ).
Example:
```
test/testHardisp data/NTUA.BLQ > onsa.tmp
paste onsa.tmp test/test_onsa_results | awk '{printf "%9.6f %9.6f %9.6f\n", $1-$4, $2-$5, $3-$6}'
```

Again, the differences (per line and column) should not exceed 1e<sup>-6</sup>

### interp  <a id="interp-cmp"></a>


## gpt3

GPT3 is a new/refined empirical set of trpospheric mapping functions, described at 
Landskron & Böhm (2018) VMF3/GPT3: refined discrete and empirical troposphere mapping functions. DOI:10.1007/s00190-017-1066-2.
It is not included in the IERS2010 standards but we include a C++ version of the 
core algorithms here, adopted from [gpt3_5_fast.m](https://vmf.geo.tuwien.ac.at/codes/gpt3_5_fast.m)
and [gpt3_1_fast.m](https://vmf.geo.tuwien.ac.at/codes/gpt3_1_fast.m. 
This function is not wrapped in the iers2010 namespace, but in a namespace called 
`dso`. Users must also include the `tropo.hpp` header file.

GPT3 needs a grid file to work, which can be downloaded from TU Wienn, at
https://vmf.geo.tuwien.ac.at/codes/

### Checking gpt3 implementation

Two dedicated test cases are provided to check the constistency of the original 
Octave/Matlab implementation provided by TU Vienna, against the one implemented 
here; these are `test/test_gpt35.cpp` and `test/test_gpt31.cpp`; results from 
two implementations agree (at least) within 1e-12 in respective units.

# How to use the library <a name="how-to"></a>

## Namespaces

- __namespace `iers2010` includes all model implementation functions.__

- namespace `iers2010::dhtide` includes details and functions only relevant to 
  the `dehanttideinel` funtion. You should probably never have to use this.

- namespace `iers2010::hisp` includes details and functions only relevant to 
  the `hardisp` funtion. You should probably never have to use this.

- namespace `iers2010::oeop` includes details and functions only relevant to 
  the `orthoeop` funtion. You should probably never have to use this.


## Linking

- static
- dynamic

# Documentation & Library API (TODO) <a name="dox"></a>

- build dox with doxygen (or link to dox)

# TODO <a name="todo"></a>

- [x] ~~test compilation against c++17 (gcc)~~
- [x] ~~the new version of dehanttidenl has a new example test case; use it!~~
- [ ] the fourth case of dehanttideinel fails with large discrepancies; wtf?

Three test programs fail:
 - testCnmtx
 - testOrthoeop
 - testRgZont2

# Bugs & Maintanance <a name="bugs"></a>
Xanthos, xanthos@mail.ntua.gr
Mitsos, danast@mail.ntua.gr

-->
