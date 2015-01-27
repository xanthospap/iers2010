gfortran FUNDARG.F -o fundarg.e
gfortran -c I_FUNDARG.F -o FUNDARG.o
gfortran -c PMSDNUT2.F -o PMSDNUT2.o
gfortran -c UTLIBR.F -o UTLIBR.o
gfortran -c FCNNUT.F -o FCNNUT.o
gfortran -c ARG2.F -o ARG2.o
gfortran -c ORTHO_EOP.F -o ORTHO_EOP.o
gfortran -c CNMTX.F -o CNMTX.o
gfortran -c RG_ZONT2.F -o RG_ZONT2.o

gfortran test.f -o test.e PMSDNUT2.o UTLIBR.o FCNNUT.o ARG2.o ORTHO_EOP.o RG_ZONT2.o CNMTX.o FUNDARG.o