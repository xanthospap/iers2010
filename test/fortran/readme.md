## Creating test data for test/iers

This folder contains scripts and source code to crete test data (some of) the 
tests residing in `test/iers`. To create the datasets, you need to:

* run the script `prepareFortranTests.sh` which will download source code files 
  published by the IERS (in FORTRAN).
* compile the test programs (linking against the IERS source code); you just 
  need to tigger `make`.
* run the individual programs, and redirect output to (assuming you are in the 
  folder `test/fortran`):
    - `./TEST_ORTHO_EOP > ../iers/data/ortho_eop.dat`

