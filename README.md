# libiers10++

[![Build Status](https://travis-ci.com/xanthospap/ggeodesy.svg?branch=master)](https://travis-ci.com/xanthospap/ggeodesy)

C++ library implementing the IERS 2010 standards.

## Introduction
This project contains a number of functions implementing models defined in
IERS Conventions (2010). The functions are available in FORTRAN from the [IERS website](https://iers-conventions.obspm.fr/conventions_material.php). Note that the
software found at this website is routinely updated.
The FORTRAN subroutines are translated to C++ with (as much as
possible) minor modifications.
__Note that the Software is now (January 2017) available at the [ftp site](ftp://maia.usno.navy.mil/conventions/2010/2010_update/) in the `software` folder per chapter.__ 

## Prerequisites

The C++ library [ggdatetime](https://github.com/xanthospap/ggdatetime) is used in the library functions to handle datetime instances when needed. Hence, you should
have [ggdatetime](https://github.com/xanthospap/ggdatetime) on your system.

Other than that, you will need a C++ compiler and (at least at this point) the `autoreconf` program which is part of the
[GNU Autotools](https://en.wikipedia.org/wiki/GNU_Autotools).

## Compilation / Installation

Source code is ISO C++17. Compilation should be trivial using any gcc version 
supporting the c++17 standard (option `-std=c++17`).

> This software is meant to be implemented on Unix-type OS's. No effort will be
> undertaken for compatibility with other OS types.

To compile the library, just follow the basic steps: (*note that the library is still at development phase so users need to configure the project before compiling*)

**If you do not need the DEBUG version** (which most probably you don't), create the `Makefile.am` templates. This means that you
should rename [Makefile.am.production](src/Makefile.am.production) and [Makefile.am.production](test/Makefile.am.production) to
`src/Makefile.am` and `test/Makefile.am` respectively.

Then run Autotools and compile:

```
autoreconf -if
./configure
make
make install
```

### Distributed Data Files

> The data file `gpt2_5.grd` is only needed for (i.e. is read by) the function `gpt2`

This repository includes also the data file [gpt2_5.grd](data/gpt2_5.grd) in the
`data` directory. This file is needed for computations when the function [gpt2](#gpt2-cmp) in invoked.
When installing the libary (aka `make install`) this file will be installed in the 
default `share` directory of the system, under the folder `iers10`. In most X systems, 
this means that you'll end up with the file: `/usr/local/share/iers10/gpt2_5.grd`.

If you want to change the data directory path, you will need to alter the respective 
Makefile template, that is [Makefile.am](data/Makefile.am).

### BLQ format files

The library includes a helper class, i.e. `iers2010::BlqIn` (declared in [blqstrm.hpp](src/blqstrm.hpp)) 
to assist reading records off from a BLQ file. Obviously, users can make use of this 
code independent of the (rest of the) library.

The file [test_blq.cpp](test/test_blq.cpp) includes a test case usage of the 
source code for reading and manipulating BLQ files; the source code is compiled to 
the executable `testBlq`. Should you want to play with it, change the first line 
of code, aka:
```
BlqIn blq("/home/xanthos/Software/iers2010/data/NTUA.BLQ");
```
to reflect a valid BLQ file (normaly such a file is distributed within the project 
under the `/data` directory.

## Status

| Chapter | (Sub)Routine | Translated | Tested | Version  | Comments |
|:--------|:-------------|:----------:|:------:|:---------|:---------|
| 4       | [GCONV2](/https://iers-conventions.obspm.fr/content/chapter4/GCONV2.F)              |<ul><li>- [ ] </li></ul>|<ul><li>- [ ] </li></ul>| | see [ngpt car2ell](https://github.com/xanthospap/ngpt/blob/master/src/car2ell.hpp)|
| 5       | [PMSDNUT2](/https://iers-conventions.obspm.fr/content/convupdt/chapter5/PMSDNUT2.F) |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 13.11.2011 | see [pmsdnut2](#pmsdnut2-cmp) Discrepancies between FORTRAN and C++ implementation <1e<sup>-6</sup> mas; this is due to the D0 decleration|
|         | [UTLIBR](/https://iers-conventions.obspm.fr/content/chapter5/UTLIBR.F)              |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 23.06.2010 | see [utlibr](#utlibr-cmp) Discrepancies between FORTRAN and C++ implementation ~1e-7 or less; this is due to the     D0 decleration|
|         | [FUNDARG](/https://iers-conventions.obspm.fr/content/chapter5/FUNDARG.F)            |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 25.02.2010 | see [fundarg](#fundarg-cmp) C++ and FORTRAN implementations give identical results,which do not agree though with the test case provided |
|         | [FCNNUT](/https://iers-conventions.obspm.fr/content/convupdt/chapter5/FCNNUT.F)     |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 19.12.2013 | |
| 7       | [DEHANTTIDEINEL](/https://iers-conventions.obspm.fr/content/convupdt/chapter7/dehanttideinel/DEHANTTIDEINEL.F) |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 19.12.2016 | see [dehanttideinel](#dehanttideinel-cmp). One test case does not produce the expected results both in the FORTRAN and the C++ implementation; that is test case 4, as privided in the source code|
|         | [HARDISP](/https://iers-conventions.obspm.fr/content/convupdt/chapter7/hardisp/HARDISP.F) |<ul><li>- [ ] </li></ul>|<ul><li>- [ ] </li></ul>| | |
|         | [ARG2](/https://iers-conventions.obspm.fr/content/convupdt/chapter7/ARG2.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 07.10.2011 | see [arg2](#arg2-cmp)|
|         | [IERS_CMP_2015](/https://iers-conventions.obspm.fr/content/convupdt/chapter7/IERS_CMP_2015.F) | <ul><li>- [ ] </li></ul>|<ul><li>- [ ] </li></ul>| | |
| 8       | [CNMTX](/https://iers-conventions.obspm.fr/content/chapter8/CNMTX.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 17.03.2010 | see [cnmtx](#cnmtx-cmp) |
|         | [ORTHO_EOP](/https://iers-conventions.obspm.fr/content/chapter8/ORTHO_EOP.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 19.05.2010 | see [ortho_eop](#ortho_eop-cmp) |
|         | [RG_ZONT2](/https://iers-conventions.obspm.fr/content/convupdt/chapter8/RG_ZONT2.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 20.12.2011 | see [rg_zont2](#rg_zont2-cmp) |
| 9       | [FCUL_A](/https://iers-conventions.obspm.fr/content/chapter9/FCUL_A.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 13.08.2009 | see [fcul_a](#fcul_a-cmp) |
|         | [FCUL_B](/https://iers-conventions.obspm.fr/content/chapter9/FCUL_B.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 14.08.2009 | see [fcul_b](#fcul_b-cmp) |
|         | [FCUL_ZD_HPA](/https://iers-conventions.obspm.fr/content/chapter9/FCUL_ZD_HPA.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 14.08.2009 | see [fcul_zd_hpa](#fcul_zd_hpa-cmp) Both FORTRAN and C++ implementations prodice different results from the test case provided|
|         | [FCUL_ZD_HPA](/https://iers-conventions.obspm.fr/content/chapter9/FCUL_ZD_HPA.F) | <ul><li>- [x] </li></ul>|<ul><li>- [ ] </li></ul>| 14.08.2009 | see [fculzd_hpa](#fculzd_hpa-cmp) |
|         | [GMF](/https://iers-conventions.obspm.fr/content/chapter9/GMF.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 12.08.2009 | see [gmf](#gmf-cmp) |
|         | [VMF1](/https://iers-conventions.obspm.fr/content/chapter9/VMF1.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 12.01.20012 | see [vmf1](#vmf1-cmp) |
|         | [VMF1_HT](/https://iers-conventions.obspm.fr/content/chapter9/VMF1_HT.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 12.01.20012 | see [vmf1_ht](#vmf1_ht-cmp) |
|         | [GPT](/https://iers-conventions.obspm.fr/content/chapter9/GPT.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 18.10.2011 | see [gpt](#gpt-cmp) |
|         | [GPT2](/https://iers-conventions.obspm.fr/content/chapter9/GPT2.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 31.05.2013 | see [gpt2](#gpt2-cmp) |


## Test Cases & FORTRAN vs C++ implementations

> The discussion below is only a description of the discrepancies between implementations
> and does **not** depict the actual precission of the models.

## Test Programs

During the compilation, some programs are compiled to test the implementations of the individual functions in the
library. These are:

- [testIers2010](test/test_iers2010.cpp) compiled to `test/testIers2010`
  This program checks the library functions using the test cases provided in the original
  FORTRAN routines. Run and check the results; in some cases, the differences (if 
  present) are expected and reason they appear is documented in the routine-specific
  blocks below.

- [testDehantTide](test/test_dehanttide.cpp) compiled to `test/testDehantTide`
  This program checks specifically the library function `dehandtideinel` using the test cases 
  provided in the original FORTRAN routine.

- [testHardisp](test/test_hardisp.cpp) compiled to `test/testHardisp`
  This program checks the library hardisp program; to run this, you are going to need
  a `BLQ` file with records for the GNSS station 'ONSA' (or you could use the file 
  [NTUA.BLQ](data/NTUA.BLQ))


These programs contain source code to run the test cases provided in the FORTRAN implementation files.

~~If needed, alternative FORTRAN implementations are provided (in the `fortran_impl` directory) for further testing.
These can be compiled using the `fortran_impl/Makefile` aka run `make` in the `fortran_impl` folder.
Note that you will need the IERS source code (you can use the script [fetch_iers_lib.sh](fortran_impl/fetch_iers_lib.sh))
and [gfortran](https://gcc.gnu.org/fortran/). Two programs are compiled, `fortran_impl/test_iers` and
`fortran_impl/test_iers_d0`; the latter uses the alternative FORTRAN implementations (**NOT** the original
IERS source code) Alternative FORTRAN implementations are provided for:~~

~~-PMSDNUT2.F named [PMSDNUT2_D0.F](fortran_impl/PMSDNUT2_D0.F) (see [pmsdnut2](#pmsdnut2-cmp))~~
~~-UTLIBR.F named [UTLIBR_D0.F](fortran_impl/UTLIBR_D0.F) (see [utlibr](#utlibr-cmp))~~

See the individual (sub)routine chapters below for the reason these files are provided.


### fundarg <a id="fundarg-cmp"></a>

The test provided for FUNDARG (in `FUNDARG.F`) fails with maximum descripancies in the order of 
1e<sup>-11</sup> radians or less. The same thing happens when i compile and run the FORTRAN implementation
(see `fortran_impl`) exactly as provided by the IERS website. The tests for this program are coded in 
[test_fundarg.cpp](src/test_fundarg.cpp) and [test_fundarg.f](fortran_impl/test_fundarg.f).


### pmsdnut2 <a id="pmsdnut2-cmp"></a>

The test case provided (in `PMSDNUT2.F`) shows discrepancies (against the C++ implementation) in the order
of 1e<sup>-6</sup> microarcseconds or less. This is because the `DOUBLE PRECISION` float arrays `PER, XS, XC, YS and YC`
are not explicitely marked with `D0` (in `PMSDNUT2.F`). If i use the implementation 
[PMSDNUT2_D0.F](fortran_impl/PMSDNUT2_D0.F) (where
the only change from `PMSDNUT2.F` is the declerations of the arrays), then the FORTRAN and C++ 
results are identical.

*Note that the Same thing happens with [utlibr](#utlibr-cmp).*

### utlibr <a id="utlibr-cmp"></a>

The test cases provided (in `UTLIBR.F`) shows discrepancies (against the C++ implementation) in the order
of 1e<sup>-6</sup> mas and mas/day or less. This is because the `DOUBLE PRECISION` float arrays `PER, DUT1S, DUT1C, DLODS and DLODC`
are not explicitely marked with `D0` (in `UTLIBR.F`). If i use the implementation 
[UTLIBR_D0.F](fortran_impl/UTLIBR_D0.F) (where
the only change from `UTLIBR.F` is the declerations of the arrays), then the FORTRAN and C++
results are identical.

Same thing happens with [pmsdnut2](#pmsdnut2-cmp).

### fcnnut <a id="fcnnut-cmp"></a>

FORTRAN (i.e. `FCNNUT.F`) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### arg2 <a id="arg2-cmp"></a>

FORTRAN (i.e. `ARG2.F`) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### dehanttideinel <a id="dehanttideinel-cmp"></a>

There are four test cases available (in `DEHANTTIDEINEL.F`) for this routine. 
For the first three, they all yield discrepancies < 1e<sup>-13</sup> meters for the C++ 
implementation, while for the FORTRAN implementation differences are zero. __However, the 
fourth test case fails misserably, both for the C++ and the FORTRAN implementation__; in this 
case, discrepancies in the order of tens of meters are found! I do not know what 
the problem is.

### cnmtx <a id="cnmtx-cmp"></a>

FORTRAN (i.e. CNMTX.F) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### ortho_eop <a id="ortho_eop-cmp"></a>

FORTRAN (i.e. ORTHO_EOP.F) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### rg_zont2 <a id="rg_zont2-cmp"></a>

The test case (provided in RG_ZONT2.F) implemented using the C++ function, gives zero discrepancies.
**However** if i compile and run the FORTRAN source for RG_ZONT2.F, the test fails with
discrepancies up to ~5.0e-10.

### fcul_a <a id="fcul_a-cmp"></a>

FORTRAN (i.e. FCUL_A.F) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### fcul_b <a id="fcul_b-cmp"></a>

FORTRAN (i.e. FCUL_B.F) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### fcul_zd_hpa <a id="fcul_zd_hpa-cmp"></a>

FOFTRAN and C++ implementations produce identical results; however, these show discrepancies up to 4e-6 with the ones provided at the test case (in the header of FZUL_ZD_HPA.F).

### gmf <a id="gmf-cmp"></a>

FORTRAN (i.e. GMF.F) and C++ implementations produce identical results; no discrepancies found
in running the test case.


### vmf1 <a id="vmf1-cmp"></a>

FORTRAN (i.e. VMF1.F) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### vmf1_ht <a id="vmf1_ht-cmp"></a>

FORTRAN (i.e. VMF1_HT.F) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### gpt <a id="gpt-cmp"></a>

FORTRAN (i.e. GPT.F) and C++ implementations produce identical results; no discrepancies found
in running the test case.

### gpt2 <a id="gpt2-cmp"></a>

FORTRAN (i.e. GPT2.F) and C++ implementations produce discrepancies up to (i.e. max value)
~1e-9; that is way smaller than the accuracy of the model especially for the values 
P (Pressure given in hPa), T (Temperature in degrees Celsius), dT (Temperature lapse rate in degrees per km), 
E (Water vapour pressure in hPa) and UNDU (Geoid undulation in meters). For the 
mapping function coefficient (AH and AW, hydrostatic and wet respectively), the discrepancies are 
smaller than 5e-14.

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
the value 1e-6.

In the official IERS2010 software, HARDISP comes as a standalone routine/program. In
this implementation, HARDISP comes both as a standalone program as well as a function 
(see `hisp::hardisp_impl` in file [hardisp_impl.cpp](src/hardisp_impl.cpp)) that can be 
used by the users in source code.

## How to use the library (TODO)

### Namespaces

- namespace `iers2010`
- namespace `iers2010::dhtide`
- namespace `iers2010::hisp`
- namespace `iers2010::oeop`

### Linking

- static
- dynamic

## Documentation & Library API (TODO)

- build dox with doxygen (or link to dox)

## TODO

- [x] ~~test compilation against c++17 (gcc)~~
- [x] ~~the new version of dehanttidenl has a new example test case; use it!~~
- [ ] the fourth case of dehanttideinel fails with large discrepancies; wtf?

## Bugs & Maintanance
Xanthos, xanthos@mail.ntua.gr
Mitsos, danast@mail.ntua.gr

