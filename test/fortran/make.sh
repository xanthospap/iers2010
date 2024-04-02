#! /usr/bin/bash

gfortran -c FUNDARG.F
gfortran -c CNMTX.F
gfortran -c PMSDNUT2.F
gfortran -c UTLIBR.F
gfortran -c ORTHO_EOP.F
gfortran -c interp_pmut1_oceans.f
gfortran -c interp_pm_gravi.f
gfortran -o pmsdnut2.out test_pmsdnut2.f PMSDNUT2.o FUNDARG.o 
gfortran -o utlibr.out test_utlibr.f UTLIBR.o FUNDARG.o 
gfortran -o ortho_eop.out test_ortho_eop.f ORTHO_EOP.o FUNDARG.o CNMTX.o
gfortran -o eop_variations.out eop_variations.f PMSDNUT2.o UTLIBR.o ORTHO_EOP.o FUNDARG.o CNMTX.o interp_pmut1_oceans.o interp_pm_gravi.o
