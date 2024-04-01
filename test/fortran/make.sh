#! /usr/bin/bash

gfortran -c FUNDARG.F
gfortran -o pmsdnut2.out PMSDNUT2.F FUNDARG.o 
gfortran -o utlibr.out UTLIBR.F FUNDARG.o 
