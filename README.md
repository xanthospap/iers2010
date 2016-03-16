# libiers10++
C++ functions and definitions implementing the IERS 2010 standards.

## Introduction
This project contains a number of functions implementing models defined in
IERS Conventions (2010). The functions are available in FORTRAN from the [IERS website](http://maia.usno.navy.mil/conv2010/software.html). Note that the
software found at this website is routinely updated.
The FORTRAN subroutines are translated to C++ (standard c++11/c++14) with (as much as
possible) minor modifications.

## Compilation / Installation
Source code is ISO C++14. Compilation should be trivial using any gcc version 
supporting the c++14 standard (option -std=c++14).

This software is meant to be implemented on Unix-type OSs. No effort will be
undertaken for compatibility with other OS types.

To compile the library, just follow the 3 basic steps
```
./configure
make
make install
```

## Status

| Chapter | (Sub)Routine | Translated | Tested | Version  | Comments |
|:--------|:-------------|:----------:|:------:|:---------|:---------|
| 4       | [GCONV2](http://maia.usno.navy.mil/conv2010/chapter4/GCONV2.F)              |<ul><li>- [ ] </li></ul>|<ul><li>- [ ] </li></ul>| | see [ngpt car2ell](https://github.com/xanthospap/ngpt/blob/master/src/car2ell.hpp)|
| 5       | [PMSDNUT2](http://maia.usno.navy.mil/conv2010/convupdt/chapter5/PMSDNUT2.F) |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 13.11.2011 | see [pmsdnut2](#pmsdnut2-cmp) |
|         | [UTLIBR](http://maia.usno.navy.mil/conv2010/chapter5/UTLIBR.F)              |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 23.06.2010 | see [utlibr](#utlibr-cmp) |
|         | [FUNDARG](http://maia.usno.navy.mil/conv2010/chapter5/FUNDARG.F)            |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 25.02.2010 | see [fundarg](#fundarg-cmp) |
|         | [FCNNUT](http://maia.usno.navy.mil/conv2010/convupdt/chapter5/FCNNUT.F)     |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 19.12.2013 | Needs updating from IERS |
| 7       | [DEHANTTIDEINEL](http://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/DEHANTTIDEINEL.F) |<ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 24.04.2015 | see [dehanttideinel](#dehanttideinel-cmp)|
|         | [HARDISP](http://maia.usno.navy.mil/conv2010/convupdt/chapter7/hardisp/HARDISP.F) |<ul><li>- [ ] </li></ul>|<ul><li>- [ ] </li></ul>| | |
|         | [ARG2](http://maia.usno.navy.mil/conv2010/convupdt/chapter7/ARG2.F) | <ul><li>- [x] </li></ul>|<ul><li>- [x] </li></ul>| 07.10.2011 | see [arg2](#arg2-cmp)|
|         | [IERS_CMP_2015](http://maia.usno.navy.mil/conv2010/convupdt/chapter7/IERS_CMP_2015.F) | <ul><li>- [ ] </li></ul>|<ul><li>- [ ] </li></ul>| | |


## Test Cases & FORTRAN vs C++ implementations

> The discussion below is only a description of the discrepancies between implementations
> and does **not** depict the actual precission of the models.

## Test Programs

During the compilation, some programs are compiled to test the implementations of the individual functions in the
library. These are:

- [testIers2010](test/test_iers2010.cpp) compiled to `test/testIers2010` and
- [testDehantTide](test/test_dehanttide.cpp) compiled to `test/testDehantTide`

These programs contain source code to run the test cases provided in the FORTRAN implementation files.

If needed, alternative FORTRAN implementations are provided (in the `fortran_impl` directory) for further testing.
These can be compiled using the `fortran_impl/Makefile` aka run `make` in the `fortran_impl` folder.
Note that you will need the IERS source code (you can use the script [fetch_iers_lib.sh](fortran_impl/fetch_iers_lib.sh))
and [gfortran](https://gcc.gnu.org/fortran/). Two programs are compiled, `fortran_impl/test_iers` and
`fortran_impl/test_iers_d0`; the latter uses the alternative FORTRAN implementations (**NOT** the original
IERS source code) Alternative FORTRAN implementations are provided for:

- PMSDNUT2.F named [PMSDNUT2_D0.F](fortran_impl/PMSDNUT2_D0.F) (see [pmsdnut2](#pmsdnut2-cmp))
- UTLIBR.F named [UTLIBR_D0.F](fortran_impl/UTLIBR_D0.F) (see [utlibr](#utlibr-cmp))

See the individual (sub)routine chapters below for the reason these files are provided.


### fundarg <a id="fundarg-cmp"></a>

The test provided for FUNDARG (in the `FUNDARG.F` header) failed with descripancies
~1e-13 or less. The same thing happens when i compile and run the FORTRAN implementation
(see `fortran_impl`). Possibly this is because of the non-use of (explicit) scientific
format in the input arguments (e.g. `T = 0.07995893223819302` and not `T = 0.07995893223819302D0`).

[MAIN.F](fortran_impl/MAIN.F) contains the test case where all `DOUBLE PRECISION` floats
are typed using the `D0` format; the results are identical with the C++ implementation.

### pmsdnut2 <a id="pmsdnut2-cmp"></a>

The test provided shows discrepancies (against the C++ implementation in the order
of 1e-7 or less. This is because the `DOUBLE PRECISION` float arrays `PER, XS, XC, YS and YC`
are not explicitely marked with `D0` (in `PMSDNUT2.F`). If i use the implementation 
[PMSDNUT2_D0.F](fortran_impl/PMSDNUT2_D0.F) (where
the only change from `PMSDNUT2.F` is the declerations of the arrays), then the FORTRAN and C++ 
results are identical.

Same thing happens with [utlibr](#utlibr-cmp).

## utlibr <a id="utlibr-cmp"></a>

The test provided shows discrepancies (against the C++ implementation in the order
of 1e-7 or less. This is because the `DOUBLE PRECISION` float arrays `PER, DUT1S, DUT1C, DLODS and DLODC`
are not explicitely marked with `D0` (in `UTLIBR.F`). If i use the implementation 
[UTLIBR_D0.F](fortran_impl/UTLIBR_D0.F) (where
the only change from `UTLIBR.F` is the declerations of the arrays), then the FORTRAN and C++
results are identical.

Same thing happens with [pmsdnut2](#pmsdnut2-cmp).

## fcnnut <a id="fcnnut-cmp"></a>

FORTRAN (i.e. FCNNUT.F) and C++ implementations produce identical results; no discrepancies found
in running the test case.

## arg2 <a id="arg2-cmp"></a>

FORTRAN (i.e. ARG2.F) and C++ implementations produce identical results; no discrepancies found
in running the test case.

## dehanttideinel <a id="dehanttideinel-cmp"></a>

The function `dehanttideinel` is tested seperately, in the source file [test_dehanttd.cpp](test/test_dehanttd.cpp)
(build as `testDehantTide` in the `test` directory). Three individual test cases are provided,
all yieding discrepancies (between the FORTRAN and C++ implementations) smaller than ~1e-17.

## Prerequisites
None. This is a standalone library. Of course a C++ compiler is assumed!

## Bugs & Maintanance
Xanthos, xanthos@mail.ntua.gr
Mitsos, danast@mail.ntua.gr
