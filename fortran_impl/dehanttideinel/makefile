## ---------------------------------------------------------------------------
## Modified makefile to build and spit test cases for the DEHANTTIDEINEL 
## routine.
##
## WARNING! This file is not the original makefile distributed by the 
## IERS2010, it is modified!
## If you see an error for the DAT/CAL2JD routines, it might be cause in the
## corresponding source files they are declared as iau_DAT and iau_CAL2JD and
## the linker cannot find them. Change the declerations in the source files
##                                                                  xanthos
## ---------------------------------------------------------------------------

#-----------------------------------------------------------------------
#
# Description:  make file for the IERS library for DEHANTTIDEINEL subroutines.
#
# Usage:
#
#    To make the library, type:
#
#      make
#
#    To delete all object files, type:
#
#      make clean
#
# This revision:  2010 October 20 
#
# Copright (C) 2008 IERS Conventions Center.  All rights reserved.
#-----------------------------------------------------------------------

# This suite of functions is only compilable by ANSI Fortran compilers -
# give the name of your preferred Fortran compiler and compilation flags
# here.

FC = gfortran 
FF = -c -O 

#----YOU SHOULDN'T HAVE TO MODIFY ANYTHING BELOW THIS LINE---------

# Name the IERS library in its source location.

IERS_LIB = libiers-dehant.a

# The list of IERS library object files.

IERS_OBS = CAL2JD.o \
           DAT.o \
           NORM8.o \
           SPROD.o \
           ST1IDIU.o \
           ST1ISEM.o \
           ST1L1.o \
           STEP2DIU.o \
           STEP2LON.o \
           ZERO_VEC8.o \
	   DEHANTTIDEINEL.o

#EXE1 = makecc_test_cases.out
#EXE2 = check.out
EXE = dehanttide.out

#-----------------------------------------------------------------------

default: $(IERS_LIB) $(EXE)

# Make and install the library.

$(IERS_LIB): $(IERS_OBS)
	ar ru $(IERS_LIB) $?

$(EXE): $(IERS_OBS) MAIN.F
	$(FC) -o $@ MAIN.F $(IERS_OBS)

# Clean up.

clean:
	rm -f $(IERS_OBS)

#-----------------------------------------------------------------------
# The list of object file dependencies.

CAL2JD.o : CAL2JD.F 
	$(FC) $(FF) -o $@ CAL2JD.F
DAT.o : DAT.F 
	$(FC) $(FF) -o $@ DAT.F
NORM8.o : NORM8.F 
	$(FC) $(FF) -o $@ NORM8.F
SPROD.o : SPROD.F 
	$(FC) $(FF) -o $@ SPROD.F
ST1IDIU.o : ST1IDIU.F 
	$(FC) $(FF) -o $@ ST1IDIU.F
ST1ISEM.o : ST1ISEM.F 
	$(FC) $(FF) -o $@ ST1ISEM.F
ST1L1.o : ST1L1.F 
	$(FC) $(FF) -o $@ ST1L1.F
STEP2DIU.o : STEP2DIU.F 
	$(FC) $(FF) -o $@ STEP2DIU.F
STEP2LON.o : STEP2LON.F 
	$(FC) $(FF) -o $@ STEP2LON.F
ZERO_VEC8.o : ZERO_VEC8.F 
	$(FC) $(FF) -o $@ ZERO_VEC8.F
DEHANTTIDEINEL.o : DEHANTTIDEINEL.F 
	$(FC) $(FF) -o $@ DEHANTTIDEINEL.F

#-----------------------------------------------------------------------
