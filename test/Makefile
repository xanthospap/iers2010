# --------------------------------------------------------------
# MAKEFILE FOR LIBIERS10++ LIBRARY
# PART 2 / MAKEFILE FOR TEST
# test/Makefile.  Generated from Makefile.in by configure.
# --------------------------------------------------------------

##  Include the Makefile configuration file
## ===================================================
include ../Makefile.conf

## Package-related info (set by autoconf)
## ===================================================
package     = libiers10++
version     = 1.0

## VPATH-related substitution variables
## ===================================================
srcdir      = .


## Compiler-related substitution variables
## ===================================================
CXX         = g++
CXXFLAGS    = -g -O2 -std=c++11
DEFS        = -DHAVE_CONFIG_H
LDFLAGS     = 
LIBS        = 

## Compiler options
## ===================================================
CXX         = g++
DEFS        = $(OPT_LEVEL) $(USE_EXTERNAL_CONSTS) $(QUICK_EXIT)                 
CXXFLAGS    = -c -Wall $(STD_C11) $(DEFS)                                       
LDFLAGS     = -Wall $(STD_C11) $(DEFS)
LIBS        =  

## Target for this Makefile
## ===================================================
EXE         = ../test.e
SRC         = test.cpp

all:

clean:
	-rm $(EXE)

check: $(EXE)

Makefile: Makefile.in ../config.status
	cd ../ && ./config.status $@

../config.status: ../configure
	cd ../ && ./config.status --recheck

$(EXE): $(SRC) FORCE
	$(CXX) $(LDFLAGS) -I../src -L.. -o $@ $< -liers10++

FORCE:
	-rm $(EXE) 2>/dev/null

.PHONY: FORCE
