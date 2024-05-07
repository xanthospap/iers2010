#! /usr/bin/bash

gfortran -c FUNDARG.F
gfortran -c CNMTX.F
gfortran -c PMSDNUT2.F
gfortran -c UTLIBR.F
gfortran -c ORTHO_EOP.F
gfortran -c RG_ZONT2.F
gfortran -c interp_pmut1_oceans.f
gfortran -c interp_pm_gravi.f

## Ocean tide EOP variations
gfortran -o ortho_eop.out test_ortho_eop.f ORTHO_EOP.o FUNDARG.o CNMTX.o
gfortran -o eop_ocean_tide.out eop_ocean_tide.f interp_pmut1_oceans.o

## EOP variations by libration
gfortran -o eop_pmgravi.out eop_libration.f interp_pm_gravi.o UTLIBR.o FUNDARG.o
gfortran -o pmsdnut2.out test_pmsdnut2.f PMSDNUT2.o FUNDARG.o 
gfortran -o utlibr.out test_utlibr.f UTLIBR.o FUNDARG.o 
gfortran -o eop_libration.out eop_libration.f UTLIBR.o interp_pm_gravi.o FUNDARG.o

## EOP variations for plotting
gfortran -o eop_variations.out eop_variations.f interp_pm_gravi.o interp_pmut1_oceans.o UTLIBR.o FUNDARG.o

## Zonal tides on Earth's rotation
gfortran -o test-rgzont2.out test_rgzont2.f RG_ZONT2.o FUNDARG.o
