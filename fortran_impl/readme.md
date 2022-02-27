# Testing IERS2010 Library Functions

This folders contains (meta)programs, written in FORTRAN, to test results and 
descripancies between the official/original/IERS-distributed implementation and 
the C++ implementation, aka this project.

To build/run the programs contained in this folder, you will need to have the 
IERS-distributed source files from [IERS](https://iers-conventions.obspm.fr/conventions_versions.php#official_target).
Currently, the most recent version is `v 1.3.0`.

Download the Fortran source files, and use the `makefile` to build the test 
programs. **Note** that the program `dehanttideinel` has a seperate build 
folder; to buld and run see the instructions in `dehanttideinel` folder.

The programs build, are actually *meta*programs; this means that the output they 
provide is used to populate respective C++ programs to check the descripancies 
between the two implementations.

For example, to check the `pmsdnut2` program:

* (*assuming you 've build the Fortran programs*) run the program `testPmsdnut2` 

* copy the output of the program in the `test/test_pmsdnut2.cpp` source file

* set the precision via the `#define PRECISION` line in the C++ source code file
(e.g. `#define PRECISION 1e-6`); compile the C++ test programs(s) using the 
command `scons` in the top-level directory.

* run the C++ program and see descripancies

The C++ test programs need [doctest](https://github.com/doctest/doctest) to be 
build.