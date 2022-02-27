# Test DEHANTTIDEINEL

This folder contains two files; they are used to compile the original 
DEHANTITIDEINEL (list of) routines and produce a number of test cases to check
against the C++ implementation within this project.

* Download the folder `dehanttideinel`, located in the `chapter7/software dir`,
from the IERS Convention Center, distributed by the IERS

* Replace the `makefile` found therein with the `makefile` in this folder

* use `make` to build the program

* run the program `makecc_test_cases.out` to produce a set of test cases;
input values and results as obtained from the original, IERS, Fortran 
iplementation are reported on STDOUT

* copy the output of the program in the test source code file 
`test/test_dehanttideinel.cpp` and compile


## Notes

**WARNING!** This file is not the original makefile distributed by the
IERS2010, it is modified!

If you see an error for the `DAT`/`CAL2JD` routines, it might be cause in the
corresponding source files they are declared as `iau_DAT` and `iau_CAL2JD` and
the linker cannot find them. Change the declerations in the source files

xanthos 02/22
