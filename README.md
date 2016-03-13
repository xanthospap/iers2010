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


| Chapter | (Sub)Routine | Translated | Tested | Comments | Version |
|:--------|:-------------|:----------:|:------:|:---------|:--------|
| 4       | [GCONV2](http://maia.usno.navy.mil/conv2010/chapter4/GCONV2.F)|- [x]|- [x]| | |
| 5       | [PMSDNUT2](http://maia.usno.navy.mil/conv2010/convupdt/chapter5/PMSDNUT2.F) | [x] | [x] | | |
|         | [UTLIBR](http://maia.usno.navy.mil/conv2010/chapter5/UTLIBR.F)  | [x] | [x] | | |
|         | [FUNDARG](http://maia.usno.navy.mil/conv2010/chapter5/FUNDARG.F) | [x] | [x] | | |
|         | [FCNNUT](http://maia.usno.navy.mil/conv2010/convupdt/chapter5/FCNNUT.F) | [x] | [x] | | |

## Prerequisites
None. This is a standalone library. Of course a C++ compiler is assumed!

## Bugs & Maintanance
Xanthos, xanthos@mail.ntua.gr
Mitsos, danast@mail.ntua.gr
