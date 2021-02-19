# libiers10++

[![Build Status](https://travis-ci.org/xanthospap/iers2010.svg?branch=master)](https://travis-ci.org/xanthospap/iers2010)

C++ library implementing the IERS 2010 standards.

## Introduction
This project contains a number of functions implementing models defined in
[IERS Conventions (2010)](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html/) 
as described in 
[IERS Technical Note No. 36](https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36.pdf;jsessionid=716AD09DB6CA42AD5E5F2F7061EDAC0A.live2?__blob=publicationFile&v=1).
The International Earth Rotation and Reference Systems Service ([IERS](https://www.iers.org/IERS/EN/Home/home_node.html))
publishes the Conventions along with relevant documents, model implementations and 
respective test cases; the latter two are available in the FORTRAN programming 
language at the [IERS website](https://iers-conventions.obspm.fr/conventions_material.php).

This repository is an effort to translate the algorithms in the C++ programming language 
 with (as much as possible) minor modifications.
Note that the software found at this website is routinely updated.

## Prerequisites

The C++ library [ggdatetime](https://github.com/xanthospap/ggdatetime) 
is used in the library functions to handle datetime instances when needed. 
Hence, you should have [ggdatetime](https://github.com/xanthospap/ggdatetime) 
on your system.

Also *note that the library is still at development phase so users need to 
configure the project before compiling*. That is, you will need the
[GNU Autotools](https://en.wikipedia.org/wiki/GNU_Autotools) package to be 
able to install the software.

## Compilation / Installation

### TL;DR

Clone, prepare build and make!

```bash
$> git clone https://github.com/xanthospap/iers2010.git && cd iers2010
$> ./install_setup.py -c production
$> autoreconf -if
$> ./configure
$> make && sudo make install
```

### Choosing C++ Standard _(optional read)_

Source code is ISO C++17 but also __compatible with C++14__. 
Compilation should be trivial using any C++ compiler
[supporting the c++17](https://en.wikipedia.org/wiki/C%2B%2B17#Compiler_support) 
or the [c++14](https://en.wikipedia.org/wiki/C%2B%2B14#Compiler_support)
standard (option `-std=c++17` or -std=c++14` in gcc and clang). By default, the 
project build files use the C++17 standard; to specify a different one, you can either 
  * change the standard flag (`-std=c++17`) in every Makefile.am file, in the directories 
    `src`, `test` and optionally `boost`, or
  * [set the flag](#install-setup-script) when invoking the `install_setup.py` script (e.g. for C++14,
    `./install_setup.py -s 14`)

Apart from C++17 and C++14 no other standard has been tested; should you want another, 
chances are you should probably adapt the source code.

### install setup script _(optional read)_

We provide a python script ([install_setup.py](install_setup.py)) to assist the 
creation/editing of the needed Makefile.am's; this way you most probably do not need to 
mess up with any Makefiles. The script has a help message (just pass the `-h` option) 
where you can see all applicable switches. The basic options are:
  * choose between a __debug__ or a __production__ build,
  * [choose a C++ standard](#chossing_c++_standard),

For most users, just running `install_setup.py -c production` should do just fine. This will 
setup a build enviroment using the default options aka the C++17 standard, a production compilation 
mode and exclude source code depending on boost.

###  Compilation

For the following, `ROOTDIR` will be the root directory of this repository,
aka the directory under which `/src`, `/test` and `/doc` folders live.

To prepare the required files for compilation (that is the `Makefile.am` in each 
of the relevant folders) you need to run the script [install_setup.py](install_setup.py). 
You can use the `-h` switch to see the help message, but in most cases the 
command `./install_setup.py -c production` will suffice.

If needed (that is you are not running the script from `ROOTDIR`) specify the 
`ROOTDIR` path via the `-d` switch.

Then run Autotools and compile:

```bash
autoreconf -if
./configure
make
sudo make install
```

## Distributed Data Files

> The data file `gpt2_5.grd` is only needed for (i.e. is read by) the function `gpt2`

This repository includes also the data file [gpt2_5.grd](data/gpt2_5.grd) in the
`data` directory. This file is needed for computations when the function [gpt2](#gpt2-cmp) in invoked.
When installing the libary (aka `make install`) this file will be installed in the 
default `share` directory of the system, under the folder `iers10`. In most *X systems, 
this means that you'll end up with the file: `/usr/local/share/iers10/gpt2_5.grd`.

If you want to change the data directory path, you will need to alter the respective 
Makefile template, that is [data/Makefile.am](data/Makefile.am).

## BLQ format files

The library includes a helper class, i.e. `iers2010::BlqIn` (declared in 
[blqstrm.hpp](src/blqstrm.hpp)) to assist reading records off from a BLQ file. 
Obviously, users can make use of this code independent of the (rest of the) 
library.

The file [test_blq.cpp](test/test_blq.cpp) includes a test case usage of the 
source code for reading and manipulating BLQ files; the source code is compiled to 
the executable `testBlq`.

## Status of the Library _(or translated routines)_

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


## Test Cases & FORTRAN vs C++ implementations

> The discussion below is only a description of the discrepancies between implementations
> and does **not** depict the actual precission of the models.

## Test Programs _(optional read)_

To compile the test programs, you need to enter the command `make check` at the 
`ROOTDIR` folder (after you have run `make`). This will build the programs 
to test the implementations of the individual functions in the library. 
They are compiled into executables in the `test` folder:

  * `testFundarg`
  
  * `testPmsdnut2`
  
  * `testUtlibr`
  
  * `testFcnnut`
  
  * `testDehanttideinel`
  
  * `testArg2`
  
  * `testCnmtx`
  
  * `testOrthoeop`
  
  * `testRgZont2`
  
  * `testFcula`
  
  * `testFculb`
  
  * `testFculZdhPa`
  
  * `testGmf`
  
  * `testVmf1`
  
  * `testVmf1Ht`
  
  * `testGpt2` (_you need to have run `make install` first for the data file [gpt2_5.grd](#Distributed-Data-Files)
    to be in place before running this script, otherwise it will fail_)

  * `testGpt`
  
  * [testHardisp](test/test_hardisp.cpp) compiled to `test/testHardisp`
  This program checks the library hardisp program; to run this, you are going to need
  a `BLQ` file with records for the GNSS station 'ONSA' (or you could use the file 
  [NTUA.BLQ](data/NTUA.BLQ)). For details see [hardisp](#hardisp-cmp).

- ~~[testIers2010](test/test_iers2010.cpp) compiled to `test/testIers2010`
  This program checks the library functions using the test cases provided in the original
  FORTRAN routines. Run and check the results; in some cases, the differences (if 
  present) are expected and reason they appear is documented in the routine-specific
  blocks below.~~

- ~~[testDehantTide](test/test_dehanttide.cpp) compiled to `test/testDehantTide`
  This program checks specifically the library function `dehandtideinel` using the test cases 
  provided in the original FORTRAN routine.~~

These programs check the library functions against the test cases provided in 
the IERS published FORTRAN source code files. **Note that in some rare occasions, the 
original FORTRAN implementation may fail (to a given accuracy) the actual test case** 
(see e.g. [fundarg](#fundarg-cmp))

~~If needed, alternative FORTRAN implementations are provided (in the `fortran_impl` directory) for further testing.
These can be compiled using the `fortran_impl/Makefile` aka run `make` in the `fortran_impl` folder.
Note that you will need the IERS source code (you can use the script [fetch_iers_lib.sh](fortran_impl/fetch_iers_lib.sh))
and [gfortran](https://gcc.gnu.org/fortran/). Two programs are compiled, `fortran_impl/test_iers` and
`fortran_impl/test_iers_d0`; the latter uses the alternative FORTRAN implementations (**NOT** the original
IERS source code) Alternative FORTRAN implementations are provided for:~~

~~-PMSDNUT2.F named [PMSDNUT2_D0.F](fortran_impl/PMSDNUT2_D0.F) (see [pmsdnut2](#pmsdnut2-cmp))~~
~~-UTLIBR.F named [UTLIBR_D0.F](fortran_impl/UTLIBR_D0.F) (see [utlibr](#utlibr-cmp))~~

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

## How to use the library

### Namespaces

- __namespace `iers2010` includes all model implementation functions.__

- namespace `iers2010::dhtide` includes details and functions only relevant to 
  the `dehanttideinel` funtion. You should probably never have to use this.

- namespace `iers2010::hisp` includes details and functions only relevant to 
  the `hardisp` funtion. You should probably never have to use this.

- namespace `iers2010::oeop` includes details and functions only relevant to 
  the `orthoeop` funtion. You should probably never have to use this.


### Linking

- static
- dynamic

## Documentation & Library API (TODO)

- build dox with doxygen (or link to dox)

## TODO

- [x] ~~test compilation against c++17 (gcc)~~
- [x] ~~the new version of dehanttidenl has a new example test case; use it!~~
- [ ] the fourth case of dehanttideinel fails with large discrepancies; wtf?

Three test programs fail:
 - testCnmtx
 - testOrthoeop
 - testRgZont2

## Bugs & Maintanance
Xanthos, xanthos@mail.ntua.gr
Mitsos, danast@mail.ntua.gr

