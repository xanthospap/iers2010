# Matlab-Created Unit Tests for GPT3/VMF3
==============================================================================

There are three main Octave/Matlab scripts here:
  * make_tests_31.m and
  * make_tests_35.m
  * make_tests_vmf3_ht.m

These scripts will actully create source code in C++ (result is printed in 
STDOUT); this source code will contain a number of tests, that is results 
obtained from running the gpt3_[15]_fast.m and C++ code to reproduce and 
check them via the libiers2010 library.

The text created by these two scripts, should be paste (in the correct place) 
within the test/test_gpt3[15].cpp files. When compilatation for tests is 
triggered (aka 'scons --make-check) these source files will be compiled.

For these scripts to work, the folllowing files are needed (all .m files 
should be placed in the same directory):
  * gpt3_1_fast.m
  * gpt3_1_fast_readGrid.m
  * gpt3_1.grd
  * gpt3_5_fast.m
  * gpt3_5_fast_readGrid.m
  * gpt3_5.grd
  * vmf3_ht.m

Original matlab/octave source code files can be retrieved from TU Vienna, 
at https://vmf.geo.tuwien.ac.at/codes/
