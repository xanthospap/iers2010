
gfortran FUNDARG.F -o fundarg.e
gfortran -c I_FUNDARG.F -o fundarg.o
gfortran PMSDNUT2.F -o pmsdnut2.e fundarg.o
