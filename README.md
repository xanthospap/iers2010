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
configure
make
make install
```

## Status

| Chapter | (Sub)Routine | Translated | Tested | Version  | Comments |
|:--------|:-------------|:----------:|:------:|:---------|:---------|
| 4       | [GCONV2](http://maia.usno.navy.mil/conv2010/chapter4/GCONV2.F)              |<ul><li>- [ ] </li></ul>|<ul><li>- [ ] </li></ul>| | |
| 5       | [PMSDNUT2](http://maia.usno.navy.mil/conv2010/convupdt/chapter5/PMSDNUT2.F) |<ul><li>- [x] </li></ul>|<ul><li>- [ ] </li></ul>| 13.11.2011 | see [pmsdnut2](#pmsdnut2-cmp) |
|         | [UTLIBR](http://maia.usno.navy.mil/conv2010/chapter5/UTLIBR.F)              |<ul><li>- [x] </li></ul>|<ul><li>- [ ] </li></ul>| 23.06.2010 | |
|         | [FUNDARG](http://maia.usno.navy.mil/conv2010/chapter5/FUNDARG.F)            |<ul><li>- [x] </li></ul>|<ul><li>- [ ] </li></ul>| 25.02.2010 | see [fundarg](#fundarg-cmp) |
|         | [FCNNUT](http://maia.usno.navy.mil/conv2010/convupdt/chapter5/FCNNUT.F)     |<ul><li>- [x] </li></ul>|<ul><li>- [ ] </li></ul>| 19.12.2013 | Needs updating from IERS |


## Test Cases & FORTRAN vs C++ implementations

### fundarg <a id="fundarg-cmp"></a>

The test provided for FUNDARG (in the FUNDARG.F header) failed with descripancies
~1e-13 or less. The same thing happens when i compile and run the FORTRAN implementation
(see `fortran_impl`). Possibly this is because of the non-use of (explicit) scientific
format in the input arguments (e.g. `T = 0.07995893223819302` and not `T = 0.07995893223819302D0`).

`fortran_impl/MAIN.F` contains the test case where all `DOUBLE PRECISION` floats
are typed using the `D0` format; the results are identical with the C++ implementation.

### pmsdnut2 <a id="pmsdnut2-cmp"></a>

The test provided shows discrepancies (against the C++ implementation in the order
of 1e-7 or less. This is because the `DOUBLE PRECISION` float arrays `PER(J), XS(J), XC(J), YS(J) and YC(J)`
are not explicitely marked with `D0`. If i use the implementation `PMSDNUT2_D0.F` (where
the only change from `PMSDNUT2.F` is the declerations of the arrays), then the FORTRAN and C++ 
results are identical.

## Prerequisites
None. This is a standalone library. Of course a C++ compiler is assumed!

## Bugs & Maintanance
Xanthos, xanthos@mail.ntua.gr
Mitsos, danast@mail.ntua.gr
