# What is this folfer for ?

Well, it hosts FORTRAN source code mainly for testing/validation purposes, and 
secondsly for ceating plots.

# Creating Validation Data

* [test_rgzont2](test/fortran/test_rgzont2.f) creates a list of input data used 
  at the test [rgzont2](test/unit_tests/rgzont2.cpp) to validate zonal tides 
  on Earth rotation.

* [eop_ocean_tide](test/fortran/eop_ocean_tide.f) creates a list of input data used 
  at the test [eop_oceantide](test/unit_tests/eop_oceantide.cpp) to validate 
  ocean tide induced variations (diurnal/subdiurnal). It is based on the 
  IERS-published [interp](https://hpiers.obspm.fr/iers/models/interp.f) software.

* [test_ortho_eop.f](test/fortran/test_ortho_eop.f) creates a list of input data used 
  at the test [eop_oceantide_orthoeop](test/unit_tests/eop_oceantide_orthoeop.cpp) 
  to validate ocean tide induced variations (diurnal/subdiurnal). It is based on the 
  [ORTHO_EOP](https://iers-conventions.obspm.fr/content/chapter8/software/ORTHO_EOP.F) 
  software, distributed by the IERS.

# Creating Plots

## EOP variations
```
./eop_variations.out > foo & ./plot_eop_variations.py foo
```
