# Testing the C++ implementation

This folder contains C++ source code for testing the C++ implementation 
(libiers2010) against the IERS-distributed FORTRAN source code.

Two kinds of programs are provided:
  * `test_$(program_name).cpp` which include source code to test the C++ 
    implementation against the test cases provided in the distributed 
    respective FORTRAN files, and

  * `test_$(program_name)_random.cpp` which include source code for random 
    tests (computed by respective FORTAN programs) to check the C++ 
    implementation against the original FORTRAN one.

The latter test programs (`test_$(program_name)_random.cpp`), can be 
re-produced using the respective FORTRAN source code, named 
`make_test_$(program_name).f` found in the `fortran_impl` directory.

## Compilation

To compile the test programs in this folder, use the `--make-check` switch 
when invoking the `scons` build command (in the ROOT folder).

## Status

For details on the discrepancies of the C++ vs the FORTRAN implementations, see 
the source code of the programs files (`test_$(program_name)_random.cpp`). They 
assert a certain precision level, else thry fail.

| Chapter | (Sub)Routine | Tested | Revision         | Comments |
|:--------|:-------------|:------:|:-----------------|:---------|
| 8       | rg_zont2     | :heavy_check_mark:  | 2011 December 20 | |
| 8       | cnmtx        | :heavy_check_mark:  | 2010 March 17 | results agree to 1e-12 |
| 8       | ortho_eop    | :heavy_check_mark:  | 2010 March 19 | |
