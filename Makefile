# --------------------------------------------------------------
# MAKEFILE FOR PROJECT IERS2010
# --------------------------------------------------------------

#
include Makefile.conf

#
all:
	(cd cpp; make all)

test:
	(cd test; g++ -Wall -std=c++11 -L../lib/ -I../inc/ test.cpp -liers2010)
